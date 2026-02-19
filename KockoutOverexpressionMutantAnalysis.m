%% Interactive FBA: Biomass & Ethanol vs Oxygen with adjustable glucose

clear; clc;

%% 1. Load model
model = readCbModel('yeast-GEM.mat');  

%% 2. Set objective
biomassRXN = 'r_2111';   % biomass reaction
model = changeObjective(model, biomassRXN);

%% 3. Define Uptake and Output Target Reactions
glucoseRxn = 'r_1714';   % glucose uptake
oxygenRxn  = 'r_1992';   % oxygen uptake
ethanolRxn = 'r_1761';   % ethanol secretion

%% Initialize Figure
fig = figure('Name','Interactive FBA: Biomass & Ethanol','Position',[100 100 900 500]);

% Axes for biomass flux
ax1 = axes('Parent',fig,'Position',[0.1 0.55 0.35 0.4]);
xlabel(ax1,'Oxygen uptake (mmol/gDW/h)');
ylabel(ax1,'Biomass flux (h^{-1})');
title(ax1,'Biomass vs Oxygen');
grid(ax1,'on');
hold(ax1,'on');

% Axes for ethanol flux
ax2 = axes('Parent',fig,'Position',[0.55 0.55 0.35 0.4]);
xlabel(ax2,'Oxygen uptake (mmol/gDW/h)');
ylabel(ax2,'Ethanol flux (mmol/gDW/h)');
title(ax2,'Ethanol vs Oxygen');
grid(ax2,'on');
hold(ax2,'on');

% Glucose slider
uicontrol('Style','text','Position',[50 90 200 20], ...
    'String','Glucose uptake (mmol/gDW/h)');
glucoseSlider = uicontrol('Style','slider','Min',0,'Max',20,'Value',10, ...
    'Position',[50 70 200 20]);

%Oxygen slider
uicontrol('Style','text','Position',[300 90 200 20], ...
    'String','Oxygen uptake (mmol/gDW/h)');
oxygenSlider = uicontrol('Style','slider','Min',0,'Max',30,'Value',20, ...
    'Position',[300 70 200 20]);

% Show Growth Rate 
growthText = uicontrol('Style','text','Position',[600 70 200 20], ...
    'String','Growth:');

% Slider Callback
glucoseSlider.Callback = @(src,event) updatePlots( ...
    model,glucoseSlider,oxygenSlider,biomassRXN, ...
    glucoseRxn,oxygenRxn,ethanolRxn,ax1,ax2,growthText);

oxygenSlider.Callback  = @(src,event) updatePlots( ...
    model,glucoseSlider,oxygenSlider,biomassRXN, ...
    glucoseRxn,oxygenRxn,ethanolRxn,ax1,ax2,growthText);

% Initialize Plot
updatePlots(model,glucoseSlider,oxygenSlider, ...
    biomassRXN,glucoseRxn,oxygenRxn,ethanolRxn,ax1,ax2,growthText);

%% 5. Solve FBA and Plot
function updatePlots(model,glSlider,oxSlider,bioRxn,glRxn,oxRxn,ethRxn,ax1,ax2,growthText)

    % Get current slider values
    glucoseUptake = -get(glSlider,'Value');  % negative = uptake
    oxygenUptake  = -get(oxSlider,'Value');  % negative = uptake

    % Overexpression or Knockout Targets
    extraRxn1 = 'r_0163';   % Alcohol Dehydrogenase 2 Knockout. Set to zero flux.
    extraRxn2 = 'r_2115';   % Alcohol Dehydrogenase 1 Overexpression. Set to minimum 10fluxmmol/gdw/hr.
    extraRxn3 = 'r_0159';   % Acohol Acetyltransferase 1 knockout. Set to zero flux.
    % % extraRxn4 = 'r_XXXX';   % More Gene Editing Targets. 
    % extraRxn5 = 'r_XXXX';   % More Gene Editing Targets. 
    % extraRxn6 = 'r_XXXX';   % More Gene Editing Targets. 
    % Solve FBA
    modelTmp = model;
    modelTmp = changeRxnBounds(modelTmp, glRxn, glucoseUptake, 'l');
    modelTmp = changeRxnBounds(modelTmp, oxRxn, oxygenUptake, 'l');

    solution = optimizeCbModel(modelTmp,'max');

    if solution.stat == 1
        currentGrowth  = solution.f;
        currentEthanol = solution.v(strcmp(modelTmp.rxns,ethRxn));
        set(growthText,'String',sprintf('Growth: %.4f h^{-1}',currentGrowth));
    else
        currentGrowth  = NaN;
        currentEthanol = NaN;
        set(growthText,'String','Growth: infeasible');
    end

    % Sweep oxygen uptake at fixed glucose
    minO2 = 0.01;
    maxO2 = get(oxSlider,'Max');
    nPoints = 50;
    oxygenRange = logspace(log10(minO2), log10(maxO2), nPoints);

    biomassVals        = zeros(size(oxygenRange));
    ethanolVals        = zeros(size(oxygenRange));
    biomassVals_extra  = zeros(size(oxygenRange));
    ethanolVals_extra  = zeros(size(oxygenRange));

    for i = 1:length(oxygenRange)

        % Wild Type Model
        tmpModel = model;
        tmpModel = changeRxnBounds(tmpModel, glRxn, glucoseUptake, 'l');
        tmpModel = changeRxnBounds(tmpModel, oxRxn, -oxygenRange(i), 'l');
        sol = optimizeCbModel(tmpModel,'max');

        if sol.stat == 1
            biomassVals(i) = sol.f;
            ethanolVals(i) = sol.v(strcmp(tmpModel.rxns,ethRxn));
        else
            biomassVals(i) = NaN;
            ethanolVals(i) = NaN;
        end

        % Knockout model
        tmpModel = model;
        tmpModel = changeRxnBounds(tmpModel, glRxn, glucoseUptake, 'l');
        tmpModel = changeRxnBounds(tmpModel, oxRxn, -oxygenRange(i), 'l');
        tmpModel = changeRxnBounds(tmpModel, extraRxn1, 0, 'b');
        tmpModel = changeRxnBounds(tmpModel, extraRxn2, 10, 'l');
        tmpModel = changeRxnBounds(tmpModel, extraRxn3, 0, 'b');
        %tmpModel = changeRxnBounds(tmpModel, extraRxn4, 0, 'b');
        %tmpModel = changeRxnBounds(tmpModel, extraRxn5, 0, 'b');
        %tmpModel = changeRxnBounds(tmpModel, extraRxn6, 0, 'b');
        sol = optimizeCbModel(tmpModel,'max');

        if sol.stat == 1
            biomassVals_extra(i) = sol.f;
            ethanolVals_extra(i) = sol.v(strcmp(tmpModel.rxns,ethRxn));
        else
            biomassVals_extra(i) = NaN;
            ethanolVals_extra(i) = NaN;
        end
    end

    % Plot biomass vs oxygen
    cla(ax1); hold(ax1,'on');
    plot(ax1, oxygenRange, biomassVals, '-o','LineWidth',2);
    plot(ax1, oxygenRange, biomassVals_extra, '--s','LineWidth',2);
    xlabel(ax1,'Oxygen uptake (mmol/gDW/h)');
    ylabel(ax1,'Biomass flux (h^{-1})');
    title(ax1,sprintf('Biomass vs Oxygen (Glucose = %.2f)', -glucoseUptake));
    legend(ax1,'Wild Type','Knockout','Location','best');
    grid(ax1,'on');

    % Plot ethanol vs oxygen
    cla(ax2); hold(ax2,'on');
    plot(ax2, oxygenRange, ethanolVals, '-o','LineWidth',2);
    plot(ax2, oxygenRange, ethanolVals_extra, '--s','LineWidth',2);
    xlabel(ax2,'Oxygen uptake (mmol/gDW/h)');
    ylabel(ax2,'Ethanol flux (mmol/gDW/h)');
    title(ax2,sprintf('Ethanol vs Oxygen (Glucose = %.2f)', -glucoseUptake));
    legend(ax2,'Wild Type','Knockout','Location','best');
    grid(ax2,'on');

end
