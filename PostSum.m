function [T] = PostSum(post,plot)
% This is for summarying the posterior distribution samples
% Argment: post --- table format
%          plot ~ 1 = table; 0 = no table
VarNames  = post.Properties.VariableNames';
StatsNames= [{'2.5%HPD'},{'97.5%HPD'},{'2.5%PCT'},{'97.5%PCT'},{'mean'},{'median'},{'mode'},{'std'},{'skewness'},{'kurtosis'}];
stats     = nan(length(VarNames),10);
for i = 1:1:length(VarNames)
    X    = table2array(post(:,i));
    HPDI = hpdi(X,95); % 95% HPDI
    l_hp = HPDI(1);
    u_h  = HPDI(2);
    l_pc = quantile(X,0.025);
    u_pc = quantile(X,0.975);
    m    = mean(X);
    md   = median(X);
    mod  = mode2(X);
    s    = std(X);
    skw  = skewness(X);
    kts  = kurtosis(X);
    stats(i,:)= [l_hp,u_h,l_pc,u_pc,m,md,mod,s,skw,kts];
end

T = array2table(stats,'RowNames',VarNames,'VariableNames',StatsNames);
%%
if plot == 1
fig = uifigure;
uit = uitable(fig,'Data',T);
end
end

