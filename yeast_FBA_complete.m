%% Yeast-GEM FBA Analysis: Wild-Type baseline 
%% Initialize cobra toolbox before running the code 

clear; clc; close all;

fprintf('===========================================\n');
fprintf('YEAST-GEM FBA ANALYSIS\n');
fprintf('===========================================\n\n');

%% 1. Setup paths and check files
model_path = 'C:\Users\Predator\yeast-GEM\model\yeast-GEM.xml'; %%Change the path to wherever yeast-GEM model is located
bigg_dict_path = 'BiGGrxnDictionary.csv';

% Check model file
if ~exist(model_path, 'file')
    error('Model file not found: %s\nUpdate model_path in this script.', model_path);
end

% Check BiGG dictionary (This is to map the Yeast-GEM reaction ID to
% complementary BiGG reaction ID for creating escher map
if ~exist(bigg_dict_path, 'file')
    warning('BiGGrxnDictionary.csv not found! Escher JSONs will use yeast-GEM IDs.');
    use_bigg = false;
else
    use_bigg = true;
    fprintf('BiGG dictionary found: %s\n', bigg_dict_path);
end

%% 2. Load model
fprintf('Loading yeast-GEM model...\n');
model = readCbModel(model_path);
fprintf('Loaded: %d reactions, %d metabolites, %d genes\n\n', ...
    length(model.rxns), length(model.mets), length(model.genes));

%% 3. Define simulation parameters
fprintf('Setting up parameters...\n');

% Sugar exchange reactions obtained from Yeast-GEM
sugars = struct();
sugars.glucose = 'r_1714';
sugars.fructose = 'r_1709';
sugars.galactose = 'r_1713';
sugars.sucrose = 'r_2056';
sugars.maltose = 'r_1911';

sugar_names = fieldnames(sugars);

% Setting up micro-aerobic conditions
oxygen_rxn = 'r_1992';
microaerobic_O2 = -2.0;

% Sugar uptake rates
sugar_uptake_rates = [5, 10, 15, 20, 25, 30];

fprintf('Sugars: %d\n', length(sugar_names));
fprintf('Uptake rates: %d\n', length(sugar_uptake_rates));
fprintf('Total simulations: %d\n\n', length(sugar_names) * length(sugar_uptake_rates));

%% 4. Initialize results
results = struct();
for i = 1:length(sugar_names)
    sugar = sugar_names{i};
    results.(sugar).uptake = [];
    results.(sugar).biomass = [];
    results.(sugar).O2_uptake = [];
    results.(sugar).CO2_production = [];
    results.(sugar).ethanol_production = [];
    results.(sugar).yield = [];
    results.(sugar).flux_distributions = {};
end

%% 5. Run FBA simulations
fprintf('Running FBA simulations...\n\n');

for i = 1:length(sugar_names)
    sugar = sugar_names{i};
    sugar_rxn = sugars.(sugar);
    
    fprintf('Analyzing %s...\n', sugar);
    
    for j = 1:length(sugar_uptake_rates)
        uptake = sugar_uptake_rates(j);
        
        % Create model copy
        tempModel = model;
        
        % Set biomass objective
        tempModel = changeObjective(tempModel, 'r_2111');
        
        % Close all sugar uptakes
        for k = 1:length(sugar_names)
            tempModel = changeRxnBounds(tempModel, sugars.(sugar_names{k}), 0, 'l');
            tempModel = changeRxnBounds(tempModel, sugars.(sugar_names{k}), 0, 'u');
        end
        
        % Set current sugar uptake
        tempModel = changeRxnBounds(tempModel, sugar_rxn, -uptake, 'l');
        tempModel = changeRxnBounds(tempModel, sugar_rxn, 0, 'u');
        
        % Set micro-aerobic oxygen
        tempModel = changeRxnBounds(tempModel, oxygen_rxn, microaerobic_O2, 'l');
        
        % Run FBA
        solution = optimizeCbModel(tempModel, 'max');
        
        if solution.stat == 1
            % Store results
            results.(sugar).uptake(j) = uptake;
            results.(sugar).biomass(j) = solution.f;
            
            O2_idx = find(strcmp(tempModel.rxns, oxygen_rxn));
            results.(sugar).O2_uptake(j) = -solution.x(O2_idx);
            
            CO2_idx = find(strcmp(tempModel.rxns, 'r_1672'));
            results.(sugar).CO2_production(j) = solution.x(CO2_idx);
            
            EtOH_idx = find(strcmp(tempModel.rxns, 'r_1761'));
            results.(sugar).ethanol_production(j) = solution.x(EtOH_idx);
            
            MW_sugar = 180;
            results.(sugar).yield(j) = solution.f / (uptake * MW_sugar / 1000);
            
            results.(sugar).flux_distributions{j} = solution.x;
            
            fprintf('  Uptake: %.1f -> Biomass: %.4f\n', uptake, solution.f);
        else
            fprintf('  Uptake: %.1f -> FAILED\n', uptake);
            results.(sugar).uptake(j) = uptake;
            results.(sugar).biomass(j) = 0;
            results.(sugar).O2_uptake(j) = 0;
            results.(sugar).CO2_production(j) = 0;
            results.(sugar).ethanol_production(j) = 0;
            results.(sugar).yield(j) = 0;
            results.(sugar).flux_distributions{j} = [];
        end
    end
    fprintf('\n');
end

%% 6. Save results as .mat file
fprintf('Saving results...\n');
save('yeast_FBA_results.mat', 'results', 'sugar_names', 'sugar_uptake_rates', 'model');
fprintf('Saved: yeast_FBA_results.mat\n\n');

%% 7. Generate plots
fprintf('Generating plots...\n');
fig_dir = fullfile(pwd, 'figures');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
    fprintf('Created: %s\n', fig_dir);
end

colors = lines(length(sugar_names));

% Plot 1: Sugar Uptake vs Biomass
figure('Position', [100, 100, 800, 600]);
hold on;
for i = 1:length(sugar_names)
    sugar = sugar_names{i};
    plot(results.(sugar).uptake, results.(sugar).biomass, ...
        '-o', 'LineWidth', 2, 'MarkerSize', 8, 'Color', colors(i,:), 'DisplayName', sugar);
end
xlabel('Sugar Uptake (mmol/gDW/h)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Biomass Production (1/h)', 'FontSize', 12, 'FontWeight', 'bold');
title('Biomass vs Sugar Uptake (Micro-aerobic)', 'FontSize', 14);
legend('Location', 'northwest'); grid on;
saveas(gcf, fullfile(fig_dir, 'sugar_uptake_vs_biomass.png'));

% Plot 2: O2 Uptake vs Biomass
figure('Position', [100, 100, 800, 600]);
hold on;
for i = 1:length(sugar_names)
    sugar = sugar_names{i};
    plot(results.(sugar).O2_uptake, results.(sugar).biomass, ...
        '-s', 'LineWidth', 2, 'MarkerSize', 8, 'Color', colors(i,:), 'DisplayName', sugar);
end
xlabel('O_2 Uptake (mmol/gDW/h)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Biomass Production (1/h)', 'FontSize', 12, 'FontWeight', 'bold');
title('Biomass vs Oxygen Uptake (Micro-aerobic)', 'FontSize', 14);
legend('Location', 'northwest'); grid on;
saveas(gcf, fullfile(fig_dir, 'oxygen_uptake_vs_biomass.png'));

% Plot 3: Yield
figure('Position', [100, 100, 800, 600]);
hold on;
for i = 1:length(sugar_names)
    sugar = sugar_names{i};
    plot(results.(sugar).uptake, results.(sugar).yield, ...
        '-d', 'LineWidth', 2, 'MarkerSize', 8, 'Color', colors(i,:), 'DisplayName', sugar);
end
xlabel('Sugar Uptake (mmol/gDW/h)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Biomass Yield (g/g)', 'FontSize', 12, 'FontWeight', 'bold');
title('Biomass Yield (Micro-aerobic)', 'FontSize', 14);
legend('Location', 'best'); grid on;
saveas(gcf, fullfile(fig_dir, 'biomass_yield.png'));

% Plot 4: Max Biomass Comparison
figure('Position', [100, 100, 800, 600]);
max_biomass = zeros(length(sugar_names), 1);
for i = 1:length(sugar_names)
    sugar = sugar_names{i};
    max_biomass(i) = max(results.(sugar).biomass);
end
bar(max_biomass, 'FaceColor', 'flat', 'CData', colors);
set(gca, 'XTickLabel', sugar_names);
ylabel('Maximum Biomass Production (1/h)', 'FontSize', 12, 'FontWeight', 'bold');
title('Maximum Biomass Comparison', 'FontSize', 14);
grid on; xtickangle(45);
saveas(gcf, fullfile(fig_dir, 'max_biomass_comparison.png'));

% Plot 5: Ethanol Production
figure('Position', [100, 100, 800, 600]);
hold on;
for i = 1:length(sugar_names)
    sugar = sugar_names{i};
    plot(results.(sugar).uptake, results.(sugar).ethanol_production, ...
        '-^', 'LineWidth', 2, 'MarkerSize', 8, 'Color', colors(i,:), 'DisplayName', sugar);
end
xlabel('Sugar Uptake (mmol/gDW/h)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Ethanol Production (mmol/gDW/h)', 'FontSize', 12, 'FontWeight', 'bold');
title('Ethanol Production (Micro-aerobic)', 'FontSize', 14);
legend('Location', 'northwest'); grid on;
saveas(gcf, fullfile(fig_dir, 'ethanol_production.png'));

% Plot 6: CO2 Production
figure('Position', [100, 100, 800, 600]);
hold on;
for i = 1:length(sugar_names)
    sugar = sugar_names{i};
    plot(results.(sugar).uptake, results.(sugar).CO2_production, ...
        '-v', 'LineWidth', 2, 'MarkerSize', 8, 'Color', colors(i,:), 'DisplayName', sugar);
end
xlabel('Sugar Uptake (mmol/gDW/h)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('CO_2 Production (mmol/gDW/h)', 'FontSize', 12, 'FontWeight', 'bold');
title('CO_2 Production (Micro-aerobic)', 'FontSize', 14);
legend('Location', 'northwest'); grid on;
saveas(gcf, fullfile(fig_dir, 'CO2_production.png'));

% Plot 7: Respiratory Quotient
figure('Position', [100, 100, 800, 600]);
hold on;
for i = 1:length(sugar_names)
    sugar = sugar_names{i};
    RQ = results.(sugar).CO2_production ./ results.(sugar).O2_uptake;
    plot(results.(sugar).uptake, RQ, ...
        '-*', 'LineWidth', 2, 'MarkerSize', 8, 'Color', colors(i,:), 'DisplayName', sugar);
end
xlabel('Sugar Uptake (mmol/gDW/h)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Respiratory Quotient (CO_2/O_2)', 'FontSize', 12, 'FontWeight', 'bold');
title('Respiratory Quotient (Micro-aerobic)', 'FontSize', 14);
legend('Location', 'best'); grid on;
yline(1, '--k', 'RQ=1', 'LineWidth', 1.5);
saveas(gcf, fullfile(fig_dir, 'respiratory_quotient.png'));

% Plot 8: Metabolic Efficiency
figure('Position', [100, 100, 800, 600]);
hold on;
for i = 1:length(sugar_names)
    sugar = sugar_names{i};
    efficiency = results.(sugar).biomass ./ results.(sugar).O2_uptake;
    plot(results.(sugar).uptake, efficiency, ...
        '-p', 'LineWidth', 2, 'MarkerSize', 8, 'Color', colors(i,:), 'DisplayName', sugar);
end
xlabel('Sugar Uptake (mmol/gDW/h)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Metabolic Efficiency (Biomass/O_2)', 'FontSize', 12, 'FontWeight', 'bold');
title('Metabolic Efficiency (Micro-aerobic)', 'FontSize', 14);
legend('Location', 'best'); grid on;
saveas(gcf, fullfile(fig_dir, 'metabolic_efficiency.png'));

close all;
fprintf('Saved: figures/\n\n');

%% 8. Export Escher JSONs for each sugar 
fprintf('Exporting Escher JSON files...\n');
if ~exist('escher_maps', 'dir')
    mkdir('escher_maps');
end

% Load BiGG dictionary
if use_bigg
    fid = fopen(bigg_dict_path, 'r');
    dict_data = textscan(fid, '%s %s', 'Delimiter', ',');
    fclose(fid);
    
    yeast_ids = dict_data{1};
    bigg_ids = dict_data{2};
    bigg_map = containers.Map(yeast_ids, bigg_ids);
    fprintf('Loaded %d BiGG mappings\n', length(yeast_ids));
else
    bigg_map = [];
end

% Export each sugar
for i = 1:length(sugar_names)
    sugar = sugar_names{i};
    flux_idx = length(results.(sugar).uptake);
    flux_dist = results.(sugar).flux_distributions{flux_idx};
    
    if isempty(flux_dist)
        continue;
    end
    
    flux_data = struct();
    unmapped_count = 0;
    
    for j = 1:length(model.rxns)
        rxn_id = model.rxns{j};
        flux_value = flux_dist(j);
        
        % Include ALL reactions with significant flux (lowered threshold)
        if abs(flux_value) < 0.01
            continue;
        end
        
        % Try to map to BiGG ID
        if ~isempty(bigg_map) && isKey(bigg_map, rxn_id)
            escher_id = bigg_map(rxn_id);
        else
            % Use yeast-GEM ID if no BiGG mapping
            escher_id = rxn_id;
            unmapped_count = unmapped_count + 1;
        end
        
        % Make valid MATLAB field name
        if ~isempty(regexp(escher_id, '^[0-9]', 'once'))
            escher_id = ['R_' escher_id];
        end
        escher_id = matlab.lang.makeValidName(escher_id);
        
        flux_data.(escher_id) = flux_value;
    end
    
    % Write JSON
    json_str = jsonencode(flux_data);
    filename = sprintf('escher_maps/%s_flux.json', sugar);
    fid = fopen(filename, 'w');
    fprintf(fid, '%s', json_str);
    fclose(fid);
    
    fprintf('  %s: %d reactions (%d unmapped to BiGG)\n', sugar, length(fieldnames(flux_data)), unmapped_count);
end

% Export detailed flux report for debugging
fprintf('\nCreating detailed flux report...\n');
for i = 1:length(sugar_names)
    sugar = sugar_names{i};
    flux_idx = length(results.(sugar).uptake);
    flux_dist = results.(sugar).flux_distributions{flux_idx};
    
    if isempty(flux_dist)
        continue;
    end
    
    % Find high flux reactions (>10 mmol/gDW/h)
    high_flux = find(abs(flux_dist) > 10);
    
    filename = sprintf('escher_maps/%s_high_flux_report.txt', sugar);
    fid = fopen(filename, 'w');
    fprintf(fid, 'HIGH FLUX REACTIONS FOR %s (>10 mmol/gDW/h)\n', upper(sugar));
    fprintf(fid, '================================================\n\n');
    
    for j = 1:length(high_flux)
        idx = high_flux(j);
        fprintf(fid, 'Flux: %.4f\n', flux_dist(idx));
        fprintf(fid, 'Yeast-GEM ID: %s\n', model.rxns{idx});
        
        % Check BiGG mapping
        if ~isempty(bigg_map) && isKey(bigg_map, model.rxns{idx})
            fprintf(fid, 'BiGG ID: %s\n', bigg_map(model.rxns{idx}));
        else
            fprintf(fid, 'BiGG ID: NOT MAPPED\n');
        end
        
        fprintf(fid, 'Name: %s\n', model.rxnNames{idx});
        fprintf(fid, 'Formula: %s\n', printRxnFormula(model, model.rxns{idx}, false));
        fprintf(fid, '\n');
    end
    
    fclose(fid);
    fprintf('  Created: %s\n', filename);
end

% Export summary
summary = struct();
summary.model = 'yeast-GEM';
summary.condition = 'micro-aerobic';
summary.oxygen_mmol = -2.0;
summary.date = datestr(now);

for i = 1:length(sugar_names)
    sugar = sugar_names{i};
    summary.sugars.(sugar).max_biomass = max(results.(sugar).biomass);
    summary.sugars.(sugar).max_yield = max(results.(sugar).yield);
end

json_str = jsonencode(summary);
fid = fopen('escher_maps/summary.json', 'w');
fprintf(fid, '%s', json_str);
fclose(fid);

fprintf('\nSaved: escher_maps/\n\n');

%% 9. Done
fprintf('===========================================\n');
fprintf('ANALYSIS COMPLETE!\n');
fprintf('===========================================\n\n');
fprintf('Results:\n');
fprintf('  - yeast_FBA_results.mat\n');
fprintf('  - figures/*.png (8 plots)\n');
fprintf('  - escher_maps/*.json (flux data)\n\n');
