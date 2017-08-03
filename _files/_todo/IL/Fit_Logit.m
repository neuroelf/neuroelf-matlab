function [hrf, fit, e, param] = Fit_Logit(tc,Run,t,mode)
% function [hrf, fit, e, param] = Fit_Logit(tc,Run,t,mode)
%
% Fits IL Model, Lindquist and Wager (2007, 2009)
%
% INPUTS:
% 
% tc    - time course
% TR    - time resolution
% Run  - expermental design (indicator)
% t     - time points to estimate; length; e.g., [1:30]
% mode  - deterministic or stochastic
%   options:
%       0 - deterministic aproach 
%       1 - simulated annealing approach
%       Please note that when using simulated annealing approach you
%       may need to perform some tuning before use.
% 
% OUTPUTS:
%
% hrf   - estimated hemodynamic response function
% fit   - estimated time course
% e     - residual time course
% param - estimated amplitude, height and width %Height, time-to-peak and Width
%
% Created by Martin Lindquist on 10/02/09
%
% A simple example for one event:
% --------------------------------
% hrf = spm_hrf(.5); hrf = hrf ./ max(hrf) - .2;
% dat = hrf + .4 * rand(size(hrf));
% Run = zeros(size(dat)); Run(1) = 1;
% [h1, fit1, e1, param] = Fit_Logit(dat,Run,(1:30)',0);
% create_figure('test'); plot(dat, 'ko-', 'LineWidth', 2); hold on; plot(fit1, 'r', 'LineWidth', 2);
%
% See Example.m for more examples with multiple events.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit the Logit model

V0 = [ 1 6 1 0.5 10 1 15];                % initial values for logit fit

% Estimate theta (logit parameters)

if (mode == 1)
    disp('Stochastic Mode');
    [theta,HH,C,P,hrf,fit,e,param] = Anneal_Logit(V0,t,tc,Run);       
    
elseif (mode == 0)
    disp('Deterministic Mode');
    [theta, hrf, fit, e, param] = Det_Logit(V0,t,tc,Run);

end

end
