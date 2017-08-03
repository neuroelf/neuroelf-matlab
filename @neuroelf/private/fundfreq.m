function [f, h, p] = fundfreq(s, opts)
% fundfreq  - find fundamental frequencies in signal s
%
% FORMAT:       [f, h, p] = fundfreq(s [, opts])
%
% Input fields:
%
%       s           Sx1 signal (double)
%       opts        optional settings
%        .freq      sampling frequency of signal (default: 22050)
%        .lowf      lower bound frequency (default: 20Hz)
%        .wfunc     window function (default: 'nutall')
%
% Output fields:
%
%       f           Sx1 frequencies (number of slices)
%       h           SxH harmonics amplitudes
%       p           SxH phases
%
% See filtwin and pwelch for more information.

% persistent data for detection
persistent wresp;
if numel(wresp) ~= 1 || ...
   ~isstruct(wresp)
    wresp = struct;
    
    % create sawtooth and sinusoidal signal in 16kHz window
    st = zeros(16000, 41);
    for c = 1:41
        for t = 1:12
            to = (-1) ^ (t + 1);
            st(:, c) = st(:, c) + (to / t) .* sin((2 * pi * t * (400 + 0.1 * (c - 1)) / 16000) .* (0:15999)');
        end
        st(:, c) = st(:, c) ./ max(abs(st(:, c)));
    end
    sw = zeros(16000, 101);
    for c = 1:101
        sw(:, c) = sin((2 * pi * (400 + 0.2 * (c - 1)) / 16000) .* (0:15999)' + pi / 4);
    end
    
    % generate window shapes
    topts = struct('average', false, 'window', 4000, 'overlap', 3000, 'units', 'dB');
    wopts = struct('average', false, 'window', 800, 'overlap', 600, 'units', 'dB');
    for c = {'blackman', 'gauss', 'hamming', 'hann', 'nutall', 'rect', 'tri', 'tukey', 'welch'}
        topts.wfunc = c{1};
        wopts.wfunc = c{1};
        wout = pwelch(st, topts);
        wresp.([c{1} '_saw']) = squeeze(mean(wout(1:1500, :, :), 2));
        wout = pwelch(sw, wopts);
        wresp.([c{1} '_sin']) = squeeze(mean(wout(1:50, :, :), 2));
    end
end

% check inputs
narg = nargin;
if narg > 0 && ...
    numel(s) == 1 && ...
    isxff(s, 'wav')
    if nargin < 2 || ...
       ~isstruct(opts) || ...
        numel(opts) ~= 1
        opts = struct;
    end
    opts.freq = s.FormatChunk.SampleRate;
    s = s.Data;
    narg = max(narg, 2);
end
if narg < 1 || ...
   ~isnumeric(s) || ...
    numel(s) < 4
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
sclass = class(s);
scfac = 1;
scoff = 0;
switch (lower(sclass))
    case {'int8'}
        scfac = 1 / 128;
    case {'int16'}
        scfac = 1 / 32768;
    case {'uint8'}
        scoff = 128;
        scfac = 1 / scoff;
    case {'uint16'}
        scoff = 32768;
        scfac = 1 / scoff;
    otherwise
        if ~strcmpi(sclass, 'single') && ...
           ~strcmpi(sclass, 'double')
            error( ...
                'neuroelf:BadArgument', ...
                'Unsupported input class.' ...
            );
        end
end
if scfac ~= 1
    if scoff ~= 0
        s = scfac .* (double(s) - scoff);
    else
        s = scfac .* double(s);
    end
else
    s = limitrangec(double(s), -1, 1, 0);
end
if size(s, 2) == 2
    s = mean(s, 2);
end
s = s(:);
siglen = numel(s);
if narg < 2 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'freq') || ...
   ~isa(opts.freq, 'double') || ...
    numel(opts.freq) ~= 1 || ...
    isinf(opts.freq) || ...
    isnan(opts.freq) || ...
    opts.freq <= 0
    opts.freq = 22050;
end
if ~isfield(opts, 'lowf') || ...
   ~isa(opts.lowf, 'double') || ...
    numel(opts.lowf) ~= 1 || ...
    isinf(opts.lowf) || ...
    isnan(opts.lowf) || ...
    opts.lowf <= 0
    opts.lowf = 20;
end
if ~isfield(opts, 'wfunc') || ...
   ~ischar(opts.wfunc) || ...
   ~any(strcmpi(opts.wfunc(:)', {'blackman', 'gauss', 'hamming', ...
        'hann', 'nutall', 'rect', 'tri', 'tukey', 'welch'}))
    opts.wfunc = 'nutall';
else
    opts.wfunc = lower(opts.wfunc(:)');
end
prelen = opts.freq / opts.lowf;
slicelen = round(2 * prelen);
slicewin = filtwin(slicelen, opts.wfunc);
firstidx = 1 - slicelen + prelen;
sindices = firstidx:prelen:(siglen+slicelen-2*prelen);
tindices = (0:(slicelen-1))';
numslices = numel(sindices);

% compute frequencies for FFT bins
fftbin = opts.freq / slicelen;
freqs = 0:fftbin:(0.5*slicelen*fftbin);
nfreqs = numel(freqs);

% initialize outputs
f = NaN(numslices, 1);
if nargout > 1
    h = zeros(numslices, 20);
    if nargout > 2
        p = zeros(numslices, 20);
    end
end

% autocorrelation maximum lags
acmaxlags = min(nfreqs - 3, round(1000 / fftbin) + 2);

% normalization factor
nfactor = (1 / norm(slicewin) .^ 2);

% filter for spectrum to find resonance peaks
[sc, sf] = tempfilter(zeros(nfreqs, 1), struct('tempdct', round(500 / fftbin) + 1));
sff = (sf' * sf) \ sf';

% iterate over slice indices
for sc = 1:numslices
    
    % get data, remove mean, and multiply with window
    swin = indexarraynb(s, sindices(sc) + tindices);
    swin = swin - mean(swin);
    swin = swin .* slicewin;
    
    % compute FFT
    swin = fft(swin);
    
    % normalize and in dB
    smag = 10 .* log10(nfactor .* abs(swin(1:nfreqs)));
    
    % no signal of interest?
    ffreq = find(smag > -35);
    if isempty(ffreq)
        continue;
    end
    mmag = max(smag);
    rmag = smag - sf * (sff * smag);
    
    % plot
    plot([smag(1:200), rmag(1:200)]);
    set(gca, 'YLim', [-50, 0]);
    pause(0.1);
    
    % and peaks within 10dB or maximum
    pfreq = find(smag >= (mmag - 10) & diff([0; rmag]) >= 0 & diff([rmag; 0]) <= 0);
    
    % compute auto-correlation for a reasonable number of lags
    pacc = sxcorr(rmag, acmaxlags);
    
    % strongest peak (after first <0 value)
    pacc(1:findfirst(pacc < 0)) = 0;
    pacc(1:(findfirst(pacc > 0) - 1)) = 0;
    macc = 1 + findfirst(diff(pacc) < 0);
    
    % test peaks around resonance peak
    ptest = pfreq(:) * (1 ./ ((macc-1.5):0.02:(macc+1.5)));
    
    % fraction with good fit
    pfit = sum(abs(ptest - round(ptest)) < 0.125) ./ size(ptest, 1);
    
    % any good fit
    if any(pfit >= 0.5)
        
        % best fit
        bfit = macc + (0.02 * (findfirst(pfit == max(pfit)) - 1)) - 1.5;
        f(sc) = bfit * fftbin;        
    end
end
