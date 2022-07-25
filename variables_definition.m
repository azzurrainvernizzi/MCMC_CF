%% MCMC CF - parameters definition
f=rng;

%% Function definition to compute the next possible step in the source reagion
d_prop=@(l_stepsize,maxstep) abs(maxstep*normcdf(l_stepsize,0,1)-maxstep/2);

%% Initialization of variables
TR=1.5;% Time of repetition for fMRI info
r_min= 0.02;% minimum radius of stimulus used for VFM
radius=10.5;% Radius of the visual field display
n_samples=100;% number of iterations for MCMC
% to test this code, we set it to #100.
how_beta=1; % how_beta=0 for bCF option A,estimate beta using OLS;
%how_beta=1 for bCF option B,put beta parameter in the MCMC loop, i.e. estimate a posterior for beta
col=1;% number of column for initialize posterior
burn_in=1;% burn_in=0 do not apply burn in; burn_in=1, apply burn in
perc_burn_in=10;% percentage of iterations to ignore to allow burn in

%% Initialize parameters
p_accept = zeros(1,n_samples); % accepted solutions
VE = zeros(1,n_samples); % variance explained
post_dist = zeros(1,n_samples); % posterior distribution
loglikelihood = zeros(1,n_samples); % loglikelihood
prior_dist = zeros(1,n_samples);% prior distribution

%% Initial values of latent variables 
l_sigma=1;
proposal_width=2;
l_beta=-5;


