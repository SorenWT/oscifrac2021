This repository contains code needed to replicate the analyses in Wainio-Theberge et al. (2021): Different roles of scale-free and oscillatory activity for neural stability and flexibility. 

#Instructions:

Preprocessing of the HCP data was done using the megconnectome software (Larson-Prior et al., 2013). You can find it here: https://github.com/Washington-University/megconnectome

For preprocessing of the CamCAN data and generation of the sensor-level time-resolved IRASA spectra, please use the code from https://github.com/SorenWT/spontevo2020.

Once you have all the preprocessed  sensor-level data, you can project it into source space using SourceEst.m. This function takes a head model and source model for each individual subject. For the HCP data, these are provided with the HCP MEG data release. For the CamCAN data, they were generated from Freesurfer outputs using the scripts get_sourcemodels.m and camcan_get_headmodel_transform.m.

Fractal and oscillatory parameters were calculated using the Dynameas toolbox, which you can find here: https://github.com/SorenWT/dynameas. An example configuration structure for dm_applymeasure.m is supplied in the file config_oscifrac_main.mat. An alternative config structure config_oscifrac_fracpower.mat defines the config structure for the supplementary analyses of total fractal power. The time-resolved outputs for the analyses of intra-subject CV can be calculated using get_irasameas_timeres.m.
Please note that many of the dependencies from Dynameas are also used by the analyses functions below - please keep Dynameas in your path when running the statistics in order to avoid missing dependencies.

Finally, the statistical analyses were carried out in Oscifrac_analyses_final.m, and figures were plotted using Oscifrac_figs_final.m. These scripts can be used to generate the analyses corresponding to figures 2, 3, 5, and 6, as well as their equivalents in the supplementary materials, and the supplementary analyses of genetics in the HCP data (requires access to the restricted genetic information from the HCP data release, not provided in this repository). The analyses and figures corresponding to figure 4 can be done using oscifrac_figure4_final.m.

One final, important note: the function bootci_swt.m (used in Oscifrac_analyses_final.m) is a modified copy of MATLAB's proprietary bootci.m code; thus, it cannot be included with this repository. The code consists only of adding a third output to the function which contains the bias correction and acceleration parameters calculated in the BCA bootstrap. 