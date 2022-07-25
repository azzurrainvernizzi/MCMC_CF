%% Monte Carlo Markov Chain (MCMC) connective field model
% The MCMC CF framework is a collection of m-files that allows 
% to estimate the spatial integration properties of the
% visual system from task or resting-state data. 
% MCMC CF frameworks in based in Matlab and works with few mrVista scripts
% (https://web.stanford.edu/group/vista/cgi-bin/wiki/index.php/MrVista).
% Please note that preprocessing and data extraction needs to be done elsewhere. 
% For example purpuses, we provide processed timeseries that were analysed using mrVista.
% Current and updated versions of MCMC CF framework will be via the following link:
% http://www.visualneuroscience.nl/cf/

% How to use it:
%% Data loading
load(['/data/MCMC_CF_data/tSeries_' subject '.mat']);

%% Run MCMC CF

[data_bayes_CF] = MCMC_CF_cluster(tSeries_data)

%% Save data
save(file,'data_bayes_CF','-mat', '-v7.3');

%% Additional notes
% We welcome any observations, suggestions that may help us improve this
% toolbox. If any particular function you need is still missing, let us know
% and we'll try to incorporate it into a next release. In case you decide
% to add or modify functions yourself, please send us the modified code.
% If you think you've found a bug, please tell us. It will help greatly
% if you can supply a  minimal-length program that exhibits the bug.

%% Citing the use of the MCMC CF framwork
% We would very much appreciate it if you would (also) cite the MCMC CF
% framework preprint when you have used it in your work. 
%     Invernizzi A., Haak K.V., Carvalho J.C., Renken R.J, Cornelissen F.W.
%     'Bayesian Connective Field Modeling: a Markov Chain Monte Carlo
%     approach', doi: https://doi.org/10.1101/2020.09.03.281162
%
%   Azzurra Invernizzi
% 	1st March 2022
%	email: a.invernizzi@med.umcg.nl/azzurra.invernizzi@icloud.com
%
%	Azzurra Invernizzi, Remco Renken, Frans Cornelissen
%	Groningen, 01-03-2022
