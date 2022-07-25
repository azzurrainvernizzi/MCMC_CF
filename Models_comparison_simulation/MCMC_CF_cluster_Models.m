function [data_bayes_CF] = MCMC_CF_cluster_Models(data_sub)
variables_definition % loading all variable needed to perform MCMC_CF

idxSource = data_sub.idxSource; % Global source voxel index (in mrVista volume)
idxTarget = data_sub.idxTarget; % % Global target voxel index (in mrVista volume)
D = data_sub.D; % matrix with distances between each pair of voxel in the source region.
tSeries_source = data_sub.tSeries_source; % time serie of source region
tSeries_target = data_sub.tSeries_target; % time serie of target region
if ~isfield(data_sub,'Model')
    Model='gauss';
else
    Model=getfield(data_sub,'Model'); %#ok<*GFLD>
end

% we randomly selected our starting position based on the source voxel
% index
centerSource = idxSource(randperm(numel(idxSource),1));% find a random position in idxSource
current_centerSource = find(idxSource == centerSource); % current distance for a random position in idxSource

%% MCMC CF model
for q =1:size(tSeries_target,2) % nr of voxels across target region
    
    % we compute the MCMC CF on one time serie of target region at the time
    Y= tSeries_target(:,q);
    Y=(Y-mean(Y))/std(Y);
    varY = var(Y); % compute variance of Y
    err = varY; %#ok<*NASGU> % set error matrix to variance of Y
    
    for k=1:n_samples % number of iterations for MCMC
        
        %% bCF option A -- estimate beta using OLS
        if how_beta==0
            
            %% Latent varible: current and proposal values
            % for each latent variable and parameter that we want estimate, we propose to 'move or jump'
            % from the current position/value, elsewhere('proposal' position/value)
            l_sigma_proposal = normrnd(l_sigma,proposal_width);
            l_stepsize_proposal = normrnd(0,1);
            maxstep=max(D(:,current_centerSource)); % maxstep size should be based on source region
            st_proposal = d_prop(l_stepsize_proposal,maxstep); % inline function defined in 'variables_definition.m'
            proposal_centerSource=centerprop(st_proposal,D,current_centerSource);
            d_current=D(:,current_centerSource); % distance between current CF center and other voxels in source region
            d_proposal=D(:,proposal_centerSource);
            
            %% Computing the MCMC predictions
            % where, xi_current/_proposal contains all CF parameters (CF size)
            % VE_current/_proposal contains the variance explains
            % post/prior/loglike_current/_proposal contains the posterior,
            % prior and loglikelihood distributions respectively
            
            [xi_current,VE_current,post_current,prior_current,loglike_current] = compute_mcmc_predictions_Model(Y,tSeries_source,d_current,l_sigma,r_min,radius,how_beta,l_beta,Model);
            [xi_proposal,VE_proposal,post_proposal,prior_proposal,loglike_proposal] = compute_mcmc_predictions_Model(Y,tSeries_source,d_proposal,l_sigma_proposal,r_min,radius,how_beta,l_beta,Model); %#ok<*ASGLU>
            
            %% Acceptance probability
            p_accept(k)=exp(post_proposal-post_current);
            % if post_proposal is larger, that the probability of jumping is >1
            % therefore, we will accept this value
            testvalue = normrnd(0,1);
            accept=normcdf(testvalue,0,1);
            accepted(k)=accept<p_accept(k);
            
            %% if accepted, we update position and values of latent variables
            if accepted(k)==1
                
                current_centerSource = proposal_centerSource;
                l_sigma = l_sigma_proposal;
                VE(k) = VE_proposal;
                post_dist(k)=post_proposal;
                loglikelihood(k)=loglike_proposal;
                prior_dist(k)=prior_proposal;
                
                posterior_latent(:,col+k-1)=[l_sigma_proposal;l_beta]; % we update only the sigma latent variable while the beta is kept constant in this case
                posterior(:,col+k-1)=cat(1,xi_current,current_centerSource); % all CF parameters and CF center values are updated
                
            else
                %% otherwise, we keep position and values of latent variables
                VE(k) = VE_current;
                post_dist(k)=post_current;
                loglikelihood(k)=loglike_current;
                prior_dist(k)=prior_current;
                
                posterior_latent(:,col+k-1)=[l_sigma;l_beta];% latent variable
                posterior(:,col+k-1)=cat(1,xi_current,current_centerSource);% all CF parameters and CF center values
                
            end % end of if
            % RR edit : two lines below to keep track of AIC and BIC across
            % iterations.
            % need to change the number of free parameters depending on the model fitted.
            BICAIC(k).BIC=BIC(loglikelihood(k),numel(Y),model_parameters);%assuming one free parameter (l-sigma)
            BICAIC(k).AIC=AIC(loglikelihood(k),model_parameters);
            % end RR edit.
            
            %% for bCF option B -- put beta parameter in the MCMC loop, i.e. estimate a posterior for beta
            % the same iter describe for how_beta==0 is computed with the only
            % exception that we now include a latent variable for beta and
            % compute a posterior for it.
        elseif how_beta==1
            
            l_beta_proposal=normrnd(l_beta, proposal_width);
            l_sigma_proposal = normrnd(l_sigma,proposal_width);
            l_stepsize_proposal = normrnd(0,1);
            maxstep=max(D(:,current_centerSource));
            st_proposal = d_prop(l_stepsize_proposal,maxstep);
            proposal_centerSource=centerprop(st_proposal,D,current_centerSource);
            d_current=D(:,current_centerSource);
            d_proposal=D(:,proposal_centerSource);
            
            [xi_current,VE_current,post_current,prior_current,loglike_current] = compute_mcmc_predictions_Model(Y,tSeries_source,d_current,l_sigma,r_min,radius,how_beta,l_beta,Model);
            [xi_proposal,VE_proposal,post_proposal,prior_proposal,loglike_proposal] = compute_mcmc_predictions_Model(Y,tSeries_source,d_proposal,l_sigma_proposal,r_min,radius,how_beta,l_beta_proposal,Model);
            
            p_accept(k)=exp(post_proposal-post_current);
            testvalue = normrnd(0,1);
            accept=normcdf(testvalue,0,1);
            accepted(k)=accept<p_accept(k);
            
            
            if accepted(k)==1
                
                current_centerSource = proposal_centerSource;
                l_sigma = l_sigma_proposal;
                l_beta = l_beta_proposal;
                VE(k) = VE_proposal;
                post_dist(k)=post_proposal;
                loglikelihood(k)=loglike_proposal;
                prior_dist(k)=prior_proposal;
                
                posterior_latent(:,col+k-1)=[l_sigma_proposal;l_beta_proposal];
                posterior(:,col+k-1)=cat(1,xi_current,current_centerSource);
                
            else
                
                VE(k) = VE_current;
                post_dist(k)=post_current;
                loglikelihood(k)=loglike_current;
                prior_dist(k)=prior_current;
                
                posterior_latent(:,col+k-1)=[l_sigma;l_beta];
                posterior(:,col+k-1)=cat(1,xi_current,current_centerSource);
                
            end
            % RR edit : two lines below to keep track of AIC and BIC across
            % iterations.
            %SEE line 88 and copy the comment
            BICAIC(k).BIC=BIC(loglikelihood(k),numel(Y),model_parameters); %#ok<*AGROW> %assuming one free parameter (l-sigma)
            BICAIC(k).AIC=AIC(loglikelihood(k),model_parameters);
            %end RR edit
        end
    end
    
    
    %% Apply burn-in for each voxel of the target region (q)
    [best_fit{q},post_dist_b{q},loglikelihood_b{q},prior_dist_b{q},posterior_latent_b{q},posterior_b{q},VE_b{q}] = compute_burn_in_cf_Model (post_dist,loglikelihood,prior_dist,posterior_latent,posterior,VE,burn_in,perc_burn_in,BICAIC);
    
end

%% Storing all MCMC CF iterations for each time series as structure
% where, best_fit contains the optimal fit for all CF parameters estimated
% --> these are the values used in all our analyses and compared with the
% standard CF model
% post_dist= posterior distribution
% loglikelihood= loglikelihood distribution
% prior_dist= prior distribution
% posterior_b/posterior_latent_b= for each iteration, we store the accepted
% CF parameter and latent variable estimatations. Please note we do this as checking for our code.

% RR edit : added the BICAIC to the output structure
data_bayes_CF = struct(...
    'best_fit',best_fit...
    ,'post_dist',post_dist_b...
    ,'loglikelihood', loglikelihood_b...
    ,'prior_dist', prior_dist_b...
    ,'posterior_b',posterior_b...
    ,'posterior_latent_b', posterior_latent_b...
    ,'BICAIC',BICAIC);
% end RR edit.
end
