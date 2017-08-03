%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example code for estimating the HRF using the Inverse-Logit Model, a
% Finte Impluse Response Model and the Canonical HRF with 2 derivatives.
% Also the code illustrates our code for detecting model misspecification. 
%
% By Martin Lindquist and Tor Wager
% Edited  10/02/09
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Load time course
%

mypath = which('ilogit');
if isempty(mypath), error('Cannot find directory with ilogit.m and other functions. Not on path?'); end
[mydir] = fileparts(mypath)

load(fullfile(mydir,'timecourse'))

tc = (tc- mean(tc))/std(tc);
len = length(tc);


%% Or: create your own

randomseednum = 30; % setting this to a value will ensure replicability of random onsets

[xBF] = spm_get_bf(struct('dt', .5, 'name', 'hrf (with time and dispersion derivatives)', 'length', 32));
clear Xtrue
for i = 1:1, xx = conv(xBF.bf(:,i), [1 1 1 1 1 1 1 1 1 1 1 1 1 1]');
    Xtrue(:, i) = xx(1:66);
end
for i = 2:3, xx = conv(xBF.bf(:,i), [1 1 1 1]');
    Xtrue(:, i) = xx(1:66);
end
hrf = Xtrue * [1 .3 .2]';
xsecs = 0:.5:32;
hrf = hrf(1:length(xsecs));
hrf = hrf ./ max(hrf);
%figure; plot(xsecs, hrf, 'k')
%hrf = hrf(1:4:end); % downsample to TR, if TR is > 0.5

rand('twister', randomseednum);
R = randperm(640); R = sort(R(1:48));
Run = zeros(640,1);
for i=1:length(R), Run(R(i)) = 1; end;
true_sig = conv(Run, hrf);
true_sig = true_sig(1:640);

tc_noise = noise_arp(640, [.7 .2]);

%% SETTINGS

TR = 0.5;
T = round(30/TR);
t = 1:T;                        % samples at which to get Logit HRF Estimate
FWHM = 4;                       % FWHM for residual scan
pval = 0.01;
df = 600;
alpha = 0.001;


%%

tc = true_sig + .2 * tc_noise;

create_figure('hrf', 3, 1); 
subplot(3,1,2); axis off;
subplot(3,1,3); axis off;
subplot(3,1,1); 
han = plot(tc,'k','LineWidth', 2);
title('Sample time course'); drawnow, pause(.5);

try
    hold on;
    hh = plot_onsets(R,'k',-3,1, 1);
    drawnow
catch
    disp('Couldn''t find function to add onset sticks to plot. Skipping.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit HRF using IL-function

% Choose mode (deterministic/stochastic)

mode = 0;   % 0 - deterministic aproach 
            % 1 - simulated annealing approach
            % Please note that when using simulated annealing approach you
            % may need to perform some tuning before use.

[h1, fit1, e1, param] = Fit_Logit(tc,Run,t,mode);
[pv sres sres_ns1] = ResidScan(e1, FWHM);
[PowLoss1] = PowerLoss(e1, fit1, (len-7) , tc, TR, Run, alpha);

hold on; han(2) = plot(fit1,'r');
legend(han, {'Data' 'IL'});

drawnow, pause(.5);

disp('Summary: IL_function');

disp('Amplitude:'); disp(param(1));
disp('Time-to-peak:'); disp(param(2));
disp('Width:'); disp(param(3));

disp('MSE:'); disp((1/(len-1)*sum(e1.^2)));
disp('Mis-modeling:'); disp(pv);
disp('Power Loss:'); disp(PowLoss1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit HRF using FIR-model

% Choose mode (FIR/sFIR)

mode = 1;   % 0 - FIR 
            % 1 - smooth FIR
            
[h2, fit2, e2, param] = Fit_sFIR(tc,TR,Run,T,mode);
[pv sres sres_ns2] = ResidScan(e2, FWHM);
[PowLoss2] = PowerLoss(e2, fit2, (len-T) , tc, TR, Run, alpha);

hold on; han(3) = plot(fit2,'g');
legend(han, {'Data' 'IL' 'sFIR'});

drawnow, pause(.5);

disp('Summary: FIR');

disp('Amplitude'); disp(param(1));
disp('Time-to-peak'); disp(param(2));
disp('Width'); disp(param(3));

disp('MSE:'); disp((1/(len-1)*sum(e2.^2)));
disp('Mis-modeling'); disp(pv);
disp('Power Loss:'); disp(PowLoss2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit HRF using Canonical HRF + 2 derivatives

p=1;
      
[h3, fit3, e3, param, info] = Fit_Canonical_HRF(tc,TR,Run,T,p);
[pv sres sres_ns3] = ResidScan(e3, FWHM);
[PowLoss3] = PowerLoss(e3, fit3, (len-p) , tc, TR, Run, alpha);

hold on; han(4) = plot(fit3,'m');

legend(han,{'Data' 'IL' 'sFIR' 'DD'})
drawnow, pause(.5);

disp('Summary: Canonical + 2 derivatives');

disp('Amplitude'); disp(param(1));
disp('Time-to-peak'); disp(param(2));
disp('Width'); disp(param(3));

disp('MSE:'); disp((1/(len-1)*sum(e3.^2)));
disp('Mis-modeling'); disp(pv);
disp('Power Loss:'); disp(PowLoss3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,3); hold on; set(gca,'FontSize', 18);
title('Estimated HRF');
xsecs = 0:.5:32;
canonhrf = spm_hrf(.5); canonhrf = canonhrf ./ max(canonhrf);

han2 = plot(xsecs, hrf, 'k', 'LineWidth', 2);
hantmp = han2;
hantmp = [hantmp plot(xsecs, canonhrf, 'b', 'LineWidth', 2)];
legend(hantmp,{'True' 'Canonical'})
drawnow; pause(2);
%delete(hantmp(2));

xsecs = xsecs(1:length(h1));

han2(2) = plot(xsecs, h1,'r');
legend(han2,{'True' 'IL'})
drawnow; pause(.5);
han2(3) = plot(xsecs, h2,'g');
legend(han2,{'True' 'IL' 'sFIR'})
drawnow; pause(.5);
han2(4) = plot(xsecs, h3,'m');
legend(han2,{'True' 'IL' 'sFIR' 'DD'})
drawnow; pause(.5);

subplot(2,2,4); hold on;set(gca,'FontSize', 18);

[s1] = Fit_sFIR(sres_ns1,TR,Run,T,0);
[s2] = Fit_sFIR(sres_ns2,TR,Run,T,0);
[s3] = Fit_sFIR(sres_ns3,TR,Run,T,0);

han4 = plot(s1(1:T),'r');
hold on; han4(2) = plot(s2(1:T),'g');
hold on; han4(3) = plot(s3(1:T),'m');
hold on; plot((1:T),zeros(T,1),'--k');
legend(han4,{'IL' 'sFIR' 'DD'})
title('Mis-modeling (HRF)');

pause(2)

%% LOOP
% ******************************
% ******************************
subplot(2,2,3); cla
subplot(2,2,4); cla


for noiseperc = 0:.3:1.2
    
tc = true_sig + noiseperc * tc_noise;

create_figure('hrf', 3, 1, 1); 
% subplot(3,1,2); axis off;
% subplot(3,1,3); axis off;

subplot(3,1,1); cla; 
title('Stabiliy across varying noise levels');
han = plot(tc,'k','LineWidth', 2);
axis tight
set(gca, 'YLim', [-15 15])
title('Sample time course'); drawnow

try
    hold on;
    hh = plot_onsets(R,'k',-3,1, 1);
    drawnow
catch
    disp('Couldn''t find function to add onset sticks to plot. Skipping.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit HRF using IL-function

% Choose mode (deterministic/stochastic)

mode = 0;   % 0 - deterministic aproach 
            % 1 - simulated annealing approach
            % Please note that when using simulated annealing approach you
            % may need to perform some tuning before use.

[h1, fit1, e1, param] = Fit_Logit(tc,Run,t,mode);
[pv sres sres_ns1] = ResidScan(e1, FWHM);
[PowLoss1] = PowerLoss(e1, fit1, (len-7) , tc, TR, Run, alpha);

hold on; han(2) = plot(fit1,'r');
legend(han, {'Data' 'IL'});

drawnow

disp('Summary: IL_function');

disp('Amplitude:'); disp(param(1));
disp('Time-to-peak:'); disp(param(2));
disp('Width:'); disp(param(3));

disp('MSE:'); disp((1/(len-1)*sum(e1.^2)));
disp('Mis-modeling:'); disp(pv);
disp('Power Loss:'); disp(PowLoss1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit HRF using FIR-model

% Choose mode (FIR/sFIR)

mode = 1;   % 0 - FIR 
            % 1 - smooth FIR
            
[h2, fit2, e2, param] = Fit_sFIR(tc,TR,Run,T,mode);
[pv sres sres_ns2] = ResidScan(e2, FWHM);
[PowLoss2] = PowerLoss(e2, fit2, (len-T) , tc, TR, Run, alpha);

hold on; han(3) = plot(fit2,'g');
legend(han, {'Data' 'IL' 'sFIR'});

drawnow

disp('Summary: FIR');

disp('Amplitude'); disp(param(1));
disp('Time-to-peak'); disp(param(2));
disp('Width'); disp(param(3));

disp('MSE:'); disp((1/(len-1)*sum(e2.^2)));
disp('Mis-modeling'); disp(pv);
disp('Power Loss:'); disp(PowLoss2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit HRF using Canonical HRF + 2 derivatives

p=1;
      
[h3, fit3, e3, param, info] = Fit_Canonical_HRF(tc,TR,Run,T,p);
[pv sres sres_ns3] = ResidScan(e3, FWHM);
[PowLoss3] = PowerLoss(e3, fit3, (len-p) , tc, TR, Run, alpha);

hold on; han(4) = plot(fit3,'m');

legend(han,{'Data' 'IL' 'sFIR' 'DD'})
drawnow

disp('Summary: Canonical + 2 derivatives');

disp('Amplitude'); disp(param(1));
disp('Time-to-peak'); disp(param(2));
disp('Width'); disp(param(3));

disp('MSE:'); disp((1/(len-1)*sum(e3.^2)));
disp('Mis-modeling'); disp(pv);
disp('Power Loss:'); disp(PowLoss3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,4); cla


subplot(2,2,3); hold on; set(gca,'FontSize', 18);
title('Estimated HRF');

xsecs = 0:.5:32;
han2 = plot(xsecs, hrf, 'k', 'LineWidth', 2);
legend(han2,{'True'})
drawnow
xsecs = xsecs(1:length(h1));

han2(2) = plot(xsecs, h1,'r');
legend(han2,{'True' 'IL'})
drawnow;
han2(3) = plot(xsecs, h2,'g');
legend(han2,{'True' 'IL' 'sFIR'})
drawnow; 
han2(4) = plot(xsecs, h3,'m');
legend(han2,{'True' 'IL' 'sFIR' 'DD'})
drawnow; 

subplot(2,2,4); hold on;set(gca,'FontSize', 18);

[s1] = Fit_sFIR(sres_ns1,TR,Run,T,0);
[s2] = Fit_sFIR(sres_ns2,TR,Run,T,0);
[s3] = Fit_sFIR(sres_ns3,TR,Run,T,0);

han4 = plot(s1(1:T),'r');
hold on; han4(2) = plot(s2(1:T),'g');
hold on; han4(3) = plot(s3(1:T),'m');
hold on; plot((1:T),zeros(T,1),'--k');
legend(han4,{'IL' 'sFIR' 'DD'})
title('Mis-modeling (HRF)');
pause(1);

end