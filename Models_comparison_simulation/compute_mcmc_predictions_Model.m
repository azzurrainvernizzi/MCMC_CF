function [ xi,varE,post_dist,prior,loglikehood] = compute_mcmc_predictions_Model(Y,tSeries_source,d,l_sigma,r_min,radius,how_beta,l_beta,Model)
%added E as an output signal of this function

%% Functions definitions based on Zeidman et al. 2018
sigma = @(l_sigma,radius,r_min) (radius-r_min).*normcdf(l_sigma,0,1)+r_min;
% RR edit : added options to model flat, DoG  to the already existing SG.
% Note that for the DoG the weight for the 2nd Gaussian is fixed; the sigma
% of the 2nd Gaussian is fixed to 8 times the sigma of the first gaussian
switch lower(Model)
    case {'dog'}
        weight =  @(d,l_sigma,radius,r_min) exp((-d.^2)/(2*sigma(l_sigma,radius, r_min).^2))...
           -0.3*exp((-d.^2)/(10*sigma(l_sigma,radius, r_min).^2));
    case {'gauss','gaus','g'}
        weight = @(d,l_sigma,radius,r_min) exp((-d.^2)/(2*sigma(l_sigma,radius, r_min).^2));
    case {'disc','d','flat','f'}
        weight = @(d,l_sigma,radius,r_min) d<=sigma(l_sigma,radius, r_min); 
    otherwise
        error('invalid model selected')
end
% end RR edit
beta=@(l_beta) exp(l_beta);

%---------------------------------------------------

%% Generate predictions (convolve CF with source signals)
% rr edit 18-Jul-18
% scale the w to create a sum of 1. This will avoid the otherwise automatic
% relation between size and beta estimate. if not scaled: when size
% increases, then beta will decrease. 
w= weight(d,l_sigma,radius,r_min);
% w= w./sum(w(:));


predictions=tSeries_source*w;
% plot(predictions);pause(0.1)

if how_beta==0; % estimate beta using classical GLM and OLS
   
    varBase = ones(size(Y,1),1); % initialize VE
    X = [predictions varBase];
    B_hat = pinv(X)*Y;
    E= Y-(X*B_hat);
    varE=var(E);
    xi=[sigma(l_sigma,radius,r_min);B_hat(1)]; % storing all CF parameters: CF size (sigma).

elseif how_beta==1; % estimate a posterior for beta
    
    pTime_series_1 =beta(l_beta)*predictions;
    pTime_series_1=detrend(pTime_series_1,'constant');
    Y_demean=detrend(Y,'constant'); 
    E=Y_demean-pTime_series_1; % Azzurra
    varE=var(E);        
    xi=[sigma(l_sigma,radius,r_min);beta(l_beta)]; % storing all CF parameters: CF size (sigma) and Beta.

end

%% Computing the likelihood, prior and posterior distributions

[muhat,sigmahat] = normfit(E);
estimated_loglikelihood=log(normpdf(E,muhat,sigmahat)); % Azzurra, July 2022

loglikehood=sum(estimated_loglikelihood); % we don't include a prior for center location because it's a const
prior_s=normpdf(l_sigma,0,1);

if how_beta==0
    prior=log(prior_s);
    post_dist=sum(estimated_loglikelihood)+prior;
    
elseif how_beta==1
    prior_b=normpdf(l_beta,-2,5);
    prior=log(prior_s)+log(prior_b);
    post_dist=sum(estimated_loglikelihood)+prior;
end

end

