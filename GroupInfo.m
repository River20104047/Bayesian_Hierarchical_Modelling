function [Jj,J,N,nj] = GroupInfo(dat)
% J = Number of groups;N = sample size; nj = sample size of group_j
%   dat is a 'table', each column refers to each group
[nmax,~]  = size(dat); % group size
nj        = nmax - sum(isnan(table2array(dat))); % sample size of each group
nj(nj==0) = NaN;
J         = sum(1 - isnan(nj)); % Num of groups
Jj        = length(nj);         % nominal num of groups
N         = nansum(nj);
end