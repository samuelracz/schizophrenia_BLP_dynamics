# Supplementary code for "Diminished Variability of Alpha and Beta Band-limited Power as a Neural Signature in Schizophrenia" by Racz et al. (2025).

This folder contains the complete analysis pipeline for the aformentioned project. The scripts reproduce the main analysis outcomes and figures published in the article. Specifically,

- 'script_01_analyze_BLP_ec/eo.m' loads *pre-processed* EEG data (from eyes-closed/eyes-open conditions, respectively) and conducts the sliding window BLP          analysis, utilizing IRASA for decomposing the power spectrum into fractal and oscillatory components. Outputs (per participant) are saved into results_BLP_ec/ as Matlab workspaces.

Additionally, 'script_00_preproc_ec/eo_mara.m' provides the means to reproduce the results from raw EEG data. Raw EEG recordings will be made available at Zenodo.org upon publication. The list of included participants and the specific position of the 30s of EEG data selected for analysis is contained in the matlab workspace miscellaneous/fnames_times_ec.mat.

Note that the code loads the pre-processed EEG data from Matlab workspaces (in data_ec_mara_avg_256Hz_mat/); however, pre-processed EEG is also provided in .edf format (in data_ec_mara_avg_256Hz/).

**When using this resource, please cite** Racz, F.S., Farkas, K., Becske, M. *et al*. Reduced temporal variability of cortical excitation/inhibition ratio in schizophrenia. *Schizophr* **11**, 20 (2025). https://doi.org/10.1038/s41537-025-00568-3

Frigyes Samuel Racz

The University of Texas at Austin

email: fsr324@austin.utexas.edu

2025.
