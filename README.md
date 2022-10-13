# FTIR-filtering
Typically, FTIR (Fourier-transform infrared) absorption bands that are suitable for further analysis with chemometric methods are selected manually. This repository provides an R script to perform this pre-processing step objectively based on the reconstruction error from PCA.

More details are available in our manuscript: \
Pia Pfeiffer, Bettina Ronai, Georg Vorlaufer, Nicole Dörr, Peter Filzmoser,
Weighted LASSO variable selection for the analysis of FTIR spectra applied to the prediction of engine oil degradation, Chemometrics and Intelligent Laboratory Systems, Volume 228, 2022, 104617, ISSN 0169-7439, https://doi.org/10.1016/j.chemolab.2022.104617. (https://www.sciencedirect.com/science/article/pii/S0169743922001289) \
Abstract: The aim of this work is to quantify the relationship between different methods of artificial oil alteration as well as engine oils collected from a passenger car using FTIR (Fourier-transform infrared) spectroscopic data and chemometric methods. We propose a comprehensive procedure for the analysis of FTIR spectra: First, a reconstruction error based pre-processing to filter non-informative variables is introduced, then simultaneous variable selection and parameter estimation using the (weighted) LASSO is performed. Eventually, post-selection inference is applied to derive confidence intervals for the selected model coefficients. The proposed pre-processing methods do not rely on manual selection of FTIR absorption bands suitable for analysis but perform filtering of non-informative variables objectively. With weighted LASSO, experts’ knowledge can be integrated with the model. This pipeline for the analysis of FTIR spectroscopic data is demonstrated and validated on a real-world dataset including series of FTIR spectra of used and artificially altered engine oils. \
Keywords: High-dimensional data analysis; LASSO regression; Variable selection; Spectroscopy; Oil condition monitoring
