clear all

variables_definition; %load the other variables

SigmaNoiseTarget = 0.1; % Add a little noise to the target time series; no need to do this for the source region (real data).
MyLSigma=-0.05; %The latent parameter used to calculate sigma.
Gain=1;% A parameter to potentially change the "gain" between the source and the target. For this demonstration not used.
how_beta=1; % 0 for bCF option A (i.e. estimate beta using OLS) and 1 for bCF option B (i.e. estimate a posterior for beta).

TrueModeloption = "DoG";
TrueModelWidth = 5;
TrueCenter = 250;

switch lower(TrueModeloption)
    case "gauss"
        M=@(d)(exp((-d.^2)/(2*TrueModelWidth.^2)));
    case "dog"
        M=@(d)(exp((-d.^2)/(2*TrueModelWidth.^2))...
          -0.3*exp((-d.^2)/(10*TrueModelWidth.^2)));
otherwise
        error('wrong model selected')
end

load('tSeries_sub1.mat');

for ii = 1:10
    
    T = tSeries_data.tSeries_source * M(tSeries_data.D(:,TrueCenter)) .* Gain;
    Ttotal = zscore(T) + SigmaNoiseTarget*randn(size(T));
    Ttotal = zscore(Ttotal);
    
    tSeries_data_syn_Gauss=tSeries_data; %initialize the structure
    tSeries_data_syn_Gauss.idxTarget=tSeries_data.idxTarget(333);
    tSeries_data_syn_Gauss.tSeries_target=Ttotal;
    tSeries_data_syn_Gauss.Model='gauss';
    [Syn_bayes_CF_Gauss] = MCMC_CF_cluster_Models(tSeries_data_syn_Gauss);
    
    tSeries_data_syn_DoG=tSeries_data; %initialize the structure
    tSeries_data_syn_DoG.idxTarget=tSeries_data.idxTarget(333);
    tSeries_data_syn_DoG.tSeries_target=Ttotal;
    tSeries_data_syn_DoG.Model='dog';
    [Syn_bayes_CF_DoG] = MCMC_CF_cluster_Models(tSeries_data_syn_DoG);
    
    weight = @(d,sigma) exp((-d.^2)/(2*sigma.^2));
    w = weight(tSeries_data.D(:,Syn_bayes_CF_Gauss.best_fit(3)),Syn_bayes_CF_Gauss.best_fit(1));
    Pgauss = tSeries_data.tSeries_source * w * Syn_bayes_CF_Gauss.best_fit(2);
    
    weight = @(d,sigma) exp((-d.^2)/(2*sigma.^2))-0.3*exp((-d.^2)/(10*sigma.^2));
    w = weight(tSeries_data.D(:,Syn_bayes_CF_DoG.best_fit(3)),Syn_bayes_CF_DoG.best_fit(1));
    Pdog = tSeries_data.tSeries_source * w * Syn_bayes_CF_DoG.best_fit(2);
    
    RSS = sum((Ttotal-Pgauss).^2); SST = sum((Ttotal).^2);
    SG.Rsq(ii) = 1 - (RSS / SST);
    SG.BIC(ii) = Syn_bayes_CF_Gauss.best_fit(5);
    
    RSS = sum((Ttotal-Pdog).^2); SST = sum((Ttotal).^2);
    DoG.Rsq(ii) = 1 - (RSS / SST);
    DoG.BIC(ii) = Syn_bayes_CF_DoG.best_fit(5);
    
end

%% Plot 
figure; 
subplot(121); boxplot([SG.Rsq; DoG.Rsq]'); ylabel('Variance Exlained'); xticklabels({'SG','DoG'})
if strcmp(TrueModeloption,'Gauss'), title('True model: SG'); else, title('True model: DoG'); end
subplot(122); boxplot([SG.BIC; DoG.BIC]'); ylabel('BIC'); xticklabels({'SG','DoG'})
if strcmp(TrueModeloption,'Gauss'), title('True model: SG'); else, title('True model: DoG'); end

