function [post_table] = MCMC(logdat,para_n,para_r,para_t,para_m,para_s2,para_a,para_b,para_cj,para_dj)
% This is for MCMC calculation for hierrachical data
% It can handle missing data and nj = 1
% --------- INPUT ----------
% logdat       ~ structured data, column = group | dat is table format
% para_n       ~ number of simulations
% para_r       ~ burn-in
% para_t       ~ thinning
% para_m       ~ prior for mean of group mean
% para_s2      ~ prior for variance of group mean
% para.a       ~ prior for shape parameter of tau2  
% para.b       ~ prior for scale parameter of tau2
% para.cj      ~ prior for shape parameter of sigmaj2  
% para.dj      ~ prior for scale parameter of sigmaj2

% ---------- OUTPUT ---------
% post_table    ~ postertior samples of variable of interest
% post_mu       ~ grand mean
% post_tau2     ~ between group variance
% post_thetaj   ~ group mean
% post_sigma2j  ~ within group variance
% post_Vj       ~ variance of group means

% Extract data dimention information (n~sample size; j~group size)
[Jj,J,~,nj]= GroupInfo(logdat);

% Empty data matrix
theta_j   = nan(para_n,Jj);
sigma2_j  = nan(para_n,Jj);
mu        = nan(para_n,1);
tau2      = nan(para_n,1);
V_j       = nan(para_n,Jj); % variance of theta_j

% Initiation paramater
% (1) Group mean of #j
yij           = table2array(logdat);  
y_bar         = nanmean(yij);
njs           = nj; % for sigma2_j correction. if nj = 1, then sigma2 = Inf
njs(njs==1)   = 0;
          
theta_j(1,:)  = y_bar + rand([1,Jj]);
sigma2_j(1,:) = nanvar(table2array(logdat));
sigma2_j(sigma2_j==0) = Inf;
mu(1)         = nanmean(theta_j(1,:));
tau2(1)       = nanvar(theta_j(1,:));
V_j(1,:)      = 1 ./ (1./tau2(1) + nj./sigma2_j(1,:));

% Gibbs sampler
for i = 2:1:para_n
    % Group mean of #i
    V_j(i,:)       = 1 ./ (1./tau2(i-1) + nj./ sigma2_j(i-1,:)); 
    thetahat_j     = (mu(i-1)/tau2(i-1) + nj./sigma2_j(i-1,:).* y_bar) ./ (1./tau2(i-1) + nj./sigma2_j(i-1,:));    
    theta_j(i,:)   = normrnd(thetahat_j,sqrt(V_j(i,:)));
    
    % Grand mean
    muhat          = (nanmean(theta_j(i,:))*J/tau2(i-1) + para_m/para_s2) / (1/para_s2 + J/tau2(i-1));
    muV            = 1 / (1/para_s2 + J/tau2(i-1));
    mu(i)          = normrnd(muhat,sqrt(muV));
    
    % Within group variance
    clear shape scale
    shape          = nj ./ 2 + para_cj;
    scale          = 1./ (nansum((yij - theta_j(i,:)).^2)./2 + para_dj);
    sigma2_j(i,:)  = 1 ./ gamrnd(shape,scale);
    
    % Between group variance
    clear shape scale
    shape          = J/2 + para_a;
    scale          = 1 / nansum((theta_j(i,:) - mu(i)).^2)/2 + para_b;
    tau2(i)        = 1 ./ gamrnd(shape,scale);
end

%% Data processing
% Burn-in
mu(1:para_r)         = [];
theta_j(1:para_r,:)  = [];
sigma2_j(1:para_r,:) = [];
tau2(1:para_r)       = [];
V_j(1:para_r,:)      = [];

% Thinning
if para_t == 0 || para_t == 1
    disp('no thinning');
    post_mu         = mu;
    post_tau2       = tau2;
    post_thetaj     = theta_j;
    post_sigma2j    = sigma2_j;
    post_Vj         = V_j;
else
    post_mu         = mu(1:para_t:end);
    post_tau2       = tau2(1:para_t:end);
    post_thetaj     = theta_j(1:para_t:end,:);
    post_sigma2j    = sigma2_j(1:para_t:end,:);
    post_Vj         = V_j(1:para_t:end,:);
 
end

%% Tablize 
theta_str = 'theta';
sigma_str = 'sigma2';
V_str     = 'V';
for k = 1:1:Jj
    thetaj{k}  = sprintf('%s_%d',theta_str,k);
    sigma2j{k} = sprintf('%s_%d',sigma_str,k);
    V{k}       = sprintf('%s_%d',V_str,k);
end
VarName = [{'mu'},{'tau2'},thetaj,sigma2j,V];

post_table = array2table([post_mu,post_tau2,post_thetaj,post_sigma2j,post_Vj],'VariableNames',VarName);    

end

