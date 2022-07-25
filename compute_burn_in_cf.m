function [best_fit,post_dist_b,loglikelihood_b,prior_dist_b,posterior_latent_b,posterior_b,VE_b] = compute_burn_in_cf (post_dist,loglikelihood,prior_dist,posterior_latent,posterior,VE,burn_in,perc_burn_in)
%---------------------------------------------------
%% Functions definition

burn=@(posterior_latent,perc_burn_in) (size(posterior_latent,2)/perc_burn_in);
beta=@(l_beta) exp(l_beta);
%---------------------------------------------------
if burn_in==1 % applying burn-in -- remove the first 1000 of elements on 10000 runs
    if size(posterior_latent,2)<perc_burn_in
        burn_in=0;
    else
        n_burn=burn(posterior_latent,perc_burn_in);
        post_dist_b=post_dist(n_burn+1:end);
        posterior_latent_b= posterior_latent(:,n_burn+1:end);
        posterior_b= posterior(:,n_burn+1:end);
        loglikelihood_b=loglikelihood(:,n_burn+1:end);
        prior_dist_b=prior_dist(:,n_burn+1:end);
        VE_b=1-VE(n_burn+1:end);
        % Best fit defined using max of loglikelihood
        m=max(loglikelihood_b);
        ind=find(loglikelihood_b==m,1,'last');
        best_fit=[posterior_b(:,ind);VE(ind)];
        
        % Best fit defined using avg of loglikelihood
%         b=mean(posterior_latent_b(1,:),2);
%         best_fit=[posterior_b(:,ind);VE(ind)];
       
    end
    
elseif burn_in==0 % do not apply burn-in 
    
    post_dist_b=post_dist;
    posterior_latent_b= posterior_latent;
    posterior_b= posterior;
    loglikelihood_b=loglikelihood;
    prior_dist_b=prior_dist;
    VE_b=1-VE;
    m=max(loglikelihood_b);
    ind=find(loglikelihood_b==m,1,'last');
    best_fit=[posterior_b(:,ind);VE(ind)];
  
end
  
end

