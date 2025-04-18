# Supplementary code for "Diminished Variability of Alpha and Beta Band-limited Power as a Neural Signature in Schizophrenia" by Racz et al. (2025).

This folder contains the complete analysis pipeline for the aformentioned project. The scripts reproduce the main analysis outcomes and figures published in the manuscript. Pipeline is as follows:

- Run analysis scripts first (_script_01_analyze_BLP_ec.m_, _script_01_analyze_BLP_eo.m_, _script_01_analyze_BLP_ec_surrogate.m_, _script_01_analyze_BLP_eo_surrogate.m_). These will generate matlab workspaces for every participant separately in folders 'results_BLP_ec', 'results_BLP_eo', 'results_surrogate_ec' and 'results_surrogate_eo', storing all neural indices involved in subsequent statistical evaluation.
- - NOTE

- 

Note that the code loads the pre-processed EEG data from Matlab workspaces (in data_ec_mara_avg_256Hz_mat/); however, pre-processed EEG is also provided in .edf format (in data_ec_mara_avg_256Hz/).

**When using this resource, please cite** Racz, F.S., Farkas, K., Becske, M. *et al*. Reduced temporal variability of cortical excitation/inhibition ratio in schizophrenia. *Schizophr* **11**, 20 (2025). https://doi.org/10.1038/s41537-025-00568-3

Frigyes Samuel Racz

The University of Texas at Austin

email: fsr324@austin.utexas.edu

2025.
