# MCMC_CF
Monte Carlo Markov Chain (MCMC) connective field model

The MCMC CF framework is a collection of m-files that allows 
to estimate the spatial integration properties of the
visual system from task or resting-state data. 

It is described in:
Invernizzi A., Haak K.V., Carvalho J.C., Renken R.J, Cornelissen F.W.
 'Bayesian Connective Field Modeling: a Markov Chain Monte Carlo
 approach', doi: https://doi.org/10.1101/2020.09.03.281162

MCMC CF frameworks is Matlab-based and works with few mrVista scripts
(https://web.stanford.edu/group/vista/cgi-bin/wiki/index.php/MrVista).

Please note that preprocessing and data extraction needs to be done prior to running this code.  For example purposes, we provide processed timeseries that were analysed using mrVista.

Current and updated versions of MCMC CF framework will be accessible via the following link:
http://www.visualneuroscience.nl/cf/


Add to your path, this folder that contains: data, all needed scripts and mrVista extensions to run the MCMC CF model. 
Scripts are fully commented to allow user to follow step-by-step the code and be able to repeat it on their own dataset. 

Please note that the data are shared as Matlab files (tSeries_sub1/tSeries_sub5) that contain a structure called ‘tSeries_data’ with all data preprocessed and ready to run the MCMC CF model.  All preprocessing steps have been run in mrVista. In each structure, we stored 5 elements:

idxSource = the global source voxel indexes (in mrVista volume)
idxTarget = the global target voxel indexes (in mrVista volume)
D = matrix contains all distances between each pair of voxel in the source region. This matrix was obtained using the matlab script named ‘ccmVoxToVoxDist.m’ attached. 
tSeries_Source = time series of source region
tSeries_Target = time series of target region

