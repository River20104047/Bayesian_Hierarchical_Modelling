%% By Zijiang Yang, 2018-08-18
% Modified by 2019-09-28 
% Modified by 2019-10-10
% Modified by 2019-12-21 % Changed to normal and inverse-gamma prior

% Prepare workspace
 clc, clear, close all
 
%% Input data
dat        = readtable('A01_Data.xlsx','Sheet','Test','Range','B1:H4');
logdat     = varfun(@log,dat);
% dataname   = 'Till_Worm_T00.mat'

% MCMC parameters
para.n    = 100100; % simulation runs -  standard: 100100
para.r    = 100;  % burn in           -  stabdard: 100
para.t    = 5;    % thinning          -  standard: 5
para.s    = 10000;% resample size     -  stabdard: 10000
eff       = (para.n - para.r) / para.t

% Prior parameters
para.m    = 0;    % prior for mean of group mean          standard: 0
para.s2   = 100000;% prior for variance of group mean     standard: 10000
para.a    = 0.0001;% prior for shape parameter of tau2     standard: 0.001
para.b    = 0.0001;% prior for scale parameter of tau2     standard: 0.001
para.cj   = 0.0001;% prior for shape parameter of sigmaj2  standard: 0.001 
para.dj   = 0.0001;% prior for scale parameter of sigmaj2  standard: 0.001

%% MCMC and IC
[post]    = MCMC(logdat,para.n,para.r,para.t,para.m,para.s2,para.a,para.b,para.cj,para.dj);
% diaplot(post);                 % Diagnostic plot 
T = PostSum(post,0);             % Posterior summary table
% [AIC,BIC,DIC] = InformationCriterion(dat,post); % Imformation Criterions

%% LOO-CV
% N    = para.s;
% [Zpost,Zp,AICcv,BICcv,DICcv] = LOOCV(dat,N,para.n,para.r,para.t,para.m,para.s2,para.a,para.b,para.cj,para.dj);

%% Save variable mat
% save(dataname);
% 
% datetime










































