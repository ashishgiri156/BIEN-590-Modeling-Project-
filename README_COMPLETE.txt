========================================
YEAST-GEM FBA ANALYSIS - COMPLETE
========================================

1. Put these files in the same folder:
   - yeast_FBA_complete.m
   - BiGGrxnDictionary.csv
   - KnockoutOverexpressionMutantAnalysis.m

2. In MATLAB,
   First run: initCobraToolbox
   After that run: yeast_FBA_complete
   Finally run: KnockoutOverexpressionMutantAnalysis

3. Done.

OUTPUTS:
--------
- yeast_FBA_results.mat (all data)
- figures/*.png (8 plots)
- escher_maps/*.json (for Escher visualization)

WHAT the code does is:
-------------
- Runs FBA for 5 sugars at 6 uptake rates
- Micro-aerobic conditions (O2 = 2 mmol/gDW/h)
- Generates all plots automatically
- Exports Escher JSONs with BiGG IDs

ESCHER VISUALIZATION:
--------------------
1. Go to: https://escher.github.io/
2. Load a map (e.g., "Central metabolism")
3. Load reaction data 
4. Upload JSON from escher_maps/
5. View fluxes on pathway

========================================
