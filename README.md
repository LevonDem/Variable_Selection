# Variable_Selection
Bilevel variable selection with overlapping groups

File descriptions
-----------------
1. **Simulation_Full_Bayesian.R**: Implements variable selection for data simulated from the statistical model where Y is *categorical*. This script focuses on within-group variable selection.
2. **SampleBasedAlgorithm.cpp**: Implements our sampling-based procedure for variable selection where Y is *continuous*. The functions within this script are called from another file which first loads and processes the data to be analyzed.
3. **Diabetes_Analysis_GroupLasso.R**: Implements the penalized regression approach using different penalty functions to a diabetes dataset.
4. **Diabetes_Analysis.R**: Implements our sampling-based procedure for variable selection on a diabetes dataset.
