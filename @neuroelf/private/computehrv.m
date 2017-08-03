function opts = computehrv(onsets, opts)
% computehrv  - compute HRV from onsets
%
% FORMAT:       hrv = computehrv(onsets [, opts])
%
% Input fields:
%
%       onsets      Ox1 or 1xO vector of onsets
%       opts        optional settings
%        .detrend   detrending of time series, either 'mean' or {'linear'}
%        .diffmeth  either of 'betweenr' or {'onr'}
%        .findons   boolean flag, onsets is a Tx1 vector (default: false)
%        .fftunits  units for FFT output, either of {'squared'} or 'dB'
%        .freq      sampling frequency of onset units (default: 1)
%        .hfreq     higher HRV frequency range (default: [0.15, 0.4])
%        .lfreq     lower HRV frequency range (default: [0.04, 0.15])
%        .method    either of 'direct' or {'pwelch'}
%        .pwwin     window size for pwelch in seconds (default: 100)
%        .relpow    compute power as ratio over total power (default: true)
%        .resfreq   resampling frequency in Hertz (default: 10)
%        .usefreq   instead of computing the PSD over the RR wave form,
%                   this computes over frequency in seconds (default: false)
%
% Output fields:
%
%       hrv         struct with fields
%        .detrend   detrending selection used
%        .diffmeth  diffmeth used, either 'b' (between) or 'o' (onr)
%        .donsets   diff of onsets (in 1Hz resolution)
%        .fftunits  FFT units used, either 'd' (dB) or 's' (squared)
%        .findons   flag whether the onsets were looked up
%        .freq      given sampling frequency of onset data
%        .hfi       indices of the high-frequency range into .praw
%        .hfp       high-frequency power estimate (contribution)
%        .hfreq     selected frequency range for high frequencies
%        .lfi       indices of the low-frequency range into .praw
%        .lfp       low-frequency power estimate (contribution)
%        .lfreq     selected frequency range for low frequencies
%        .method    FFT method used, either 'd' (direct) or 'p' (pwelch)
%        .onsets    onsets (0-based, in 1Hz sampling)
%        .onsetx    estimate of positions where the BPM is given (0-based)
%        .praw      raw power spectrum (normalized FFT/pwelch output)
%        .prawf     raw power spectrum frequencies (e.g. for plot)
%        .pwo       pwelch options used (struct)
%        .pwwin     pwelch window used (only relevant if .method == 'p')
%        .relpow    flag whether .hfp/.lfp are relative to total power
%        .resfreq   resampling frequency used
%        .rmssd     RMSSD measure (alternative measure for variability)
%        .rri       RRi plot data
%        .rrifreq   inverse RRi (frequency) plot data
%        .totalp    sum of power contribution (1 if .relpow == true)
%        .totalpr   total raw power (amplitude)
%        .usefreq   flag whether computation used .rrifreq rather than .rri

% Version:  v1.0
% Build:    14091910
% Date:     Sep-19 2014, 10:43 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2014, Jochen Weber
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%     * Neither the name of Columbia University nor the
%       names of its contributors may be used to endorse or promote products
%       derived from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% argument check
if nargin < 1 || ...
   ~isa(onsets, 'double') || ...
    numel(onsets) < 3 || ...
    numel(onsets) ~= max(size(onsets)) || ...
    any(isinf(onsets) | isnan(onsets))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin < 2 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'detrend') || ...
   ~ischar(opts.detrend) || ...
    isempty(opts.detrend) || ...
   ~any(strcmpi(opts.detrend(:)', {'linear', 'mean'}))
    opts.detrend = 'linear';
else
    opts.detrend = lower(opts.detrend(:)');
end
if ~isfield(opts, 'diffmeth') || ...
   ~ischar(opts.diffmeth) || ...
    isempty(opts.diffmeth) || ...
    lower(opts.diffmeth(1)) ~= 'b'
    opts.diffmeth = 'o';
else
    opts.diffmeth = 'b';
end
if ~isfield(opts, 'fftunits') || ...
   ~ischar(opts.fftunits) || ...
   ~any(strcmpi(opts.fftunits(:)', {'d', 'db', 's', 'sq', 'squared'}))
    opts.fftunits = 's';
else
    opts.fftunits = lower(opts.units(1));
end
if ~isfield(opts, 'findons') || ...
   ~islogical(opts.findons) || ...
    numel(opts.findons) ~= 1
    opts.findons = (numel(onsets) > (numel(unique(onsets)) .^ 2));
end
hasfreq = true;
if ~isfield(opts, 'freq') || ...
   ~isa(opts.freq, 'double') || ...
    numel(opts.freq) ~= 1 || ...
    isinf(opts.freq) || ...
    isnan(opts.freq) || ...
    opts.freq <= 0
    opts.freq = 1;
    hasfreq = false;
end
if ~isfield(opts, 'hfreq') || ...
   ~isa(opts.hfreq, 'double') || ...
    numel(opts.hfreq) ~= 2 || ...
    any(isinf(opts.hfreq) | isnan(opts.hfreq) | opts.hfreq < 0) || ...
    opts.hfreq(1) >= opts.hfreq(2)
    opts.hfreq = [0.15, 0.4];
end
if ~isfield(opts, 'lfreq') || ...
   ~isa(opts.lfreq, 'double') || ...
    numel(opts.lfreq) ~= 2 || ...
    any(isinf(opts.lfreq) | isnan(opts.lfreq) | opts.lfreq < 0) || ...
    opts.lfreq(1) >= opts.lfreq(2)
    opts.lfreq = [0.04, 0.15];
end
if ~isfield(opts, 'method') || ...
   ~ischar(opts.method) || ...
    isempty(opts.method) || ...
    lower(opts.method(1)) ~= 'p'
    opts.method = 'd';
else
    opts.method = 'p';
end
if ~isfield(opts, 'pwwin') || ...
   ~isa(opts.pwwin, 'double') || ...
    numel(opts.pwwin) ~= 1 || ...
    isinf(opts.pwwin) || ...
    isnan(opts.pwwin) || ...
    opts.pwwin < (2 / min(opts.lfreq(1), opts.hfreq(1)))
    opts.pwwin = max(100, 2 / min(opts.lfreq(1), opts.hfreq(1)));
end
if ~isfield(opts, 'relpow') || ...
   ~islogical(opts.relpow) || ...
    numel(opts.relpow) ~= 1
    opts.relpow = true;
end
if ~isfield(opts, 'resfreq') || ...
   ~isa(opts.resfreq, 'double') || ...
    numel(opts.resfreq) ~= 1 || ...
    isinf(opts.resfreq) || ...
    isnan(opts.resfreq) || ...
    opts.resfreq <= 0
    opts.resfreq = 10;
end
if ~isfield(opts, 'usefreq') || ...
   ~islogical(opts.usefreq) || ...
    numel(opts.usefreq) ~= 1
    opts.usefreq = false;
end

% find onsets in signal?
if opts.findons

    % only if frequency given!
    if ~hasfreq
        error( ...
            'neuroelf:MissingOption', ...
            'For finding onsets, the .freq option must be set.' ...
        );
    end

    % locate onsets
    onsets = find(diff(onsets > mean(onsets)));

    % bail out on too short/long onset intervals
    donsets = diff(onsets);
    if any(donsets < (0.2 * opts.freq) | donsets > (3 * opts.freq))
        error( ...
            'neuroelf:BadArgument', ...
            'Unlikely onsets detected (distance too short/long).' ...
        );
    end

% otherwise
else

    % make sure onsets are in good order
    onsets = sort(onsets(:));
end

% subtract first value
onsets = onsets - onsets(1);

% convert to 1Hz if required
if opts.freq ~= 1
    onsets = (1 / opts.freq) .* onsets;
end

% what kind of diff measure
if opts.diffmeth == 'b'
    donsets = [0; 0.5 .* diff(onsets(1:end-1) + onsets(2:end))];
    donsets(1) = 2 * donsets(2) - donsets(3);
else
    donsets = diff(onsets);
end

% compute the positions we know the BPM for
onsetx = 0.5 * (onsets(1:end-1) + onsets(2:end));

% resample using regular interp
rri = interp1(onsetx, donsets, onsetx(1):(1 / opts.resfreq):onsetx(end), 'spline');
rri = rri(:);

% compute 1 / RRi
rrifreq = 1 ./ rri;

% what method to use FFT with
if opts.method == 'd'
    pwo = struct( ...
        'detrend',  opts.detrend, ...
        'nfft',     numel(rri), ...
        'overlap',  round(0.5 * numel(rri)), ...
        'units',    opts.fftunits, ...
        'window',   numel(rri));
else
    wsize = round(opts.pwwin * opts.resfreq);
    osize = round(0.5 * wsize);
    pwo = struct( ...
        'detrend',  opts.detrend, ...
        'nfft',     wsize, ...
        'overlap',  osize, ...
        'units',    opts.fftunits, ...
        'window',   wsize);
end

% compute periodogram
if opts.usefreq
    [p, w] = pwelch(rrifreq, pwo);
else
    [p, w] = pwelch(rri, pwo);
end
w = (0.5 * opts.resfreq / pi) .* w;

% compute relative power contribution
totalpr = sum(p);
lfi = find(w > opts.lfreq(1) & w <= opts.lfreq(2));
hfi = find(w > opts.hfreq(1) & w <= opts.hfreq(2));
lfp = sum(p(lfi));
hfp = sum(p(hfi));
if opts.relpow
    lfp = lfp / totalpr;
    hfp = hfp / totalpr;
    totalp = 1;
else
    totalp = totalpr;
end

% compute RMSSD
rmssd = sqrt(mean(diff(donsets) .^ 2));

% reuse opts as output!
opts.bpm = 60 .* rrifreq;
opts.bpmmean = mean(opts.bpm);
opts.donsetx = donsets;
opts.hfi = hfi;
opts.hfp = hfp;
opts.lfi = lfi;
opts.lfp = lfp;
opts.onsets = onsets;
opts.onsetx = onsetx;
opts.rmssd = rmssd;
opts.rri = rri;
opts.rrifreq = rrifreq;
opts.praw = p;
opts.prawf = w;
opts.pwo = pwo;
opts.totalp = totalp;
opts.totalpr = totalpr;
