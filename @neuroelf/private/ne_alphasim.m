function varargout = ne_alphasim(varargin)
% ne_alphasim  - run alphasim and display results in console
%
% FORMAT:       ne_alphasim(SRC, EVT, [ remote, niter, corrp, zshift, conj])
%
% Input fields:
%
%       SRC, EVT    Matlab handle callback inputs (discarded)
%       remote      for calling without interaction, set to 'remote'
%       niter       1x1 double number of iterations (set remote to 'remote')
%       corrp       1x1 double corrected p-value (set remote to 'remote')
%       zshift      1x1 double desired z-shift (set remote to 'remote')
%       conj        1x1 double number of conjunctions (remote!)
%
% Example:
%
%   ne_alphasim(0, 0, 'remote', 2500, 0.05, 0, 2);
%
% Notes: the smoothness estimate, etc. are taken from the map that is
%        selected. Only works on the first selected map!

% Version:  v1.1
% Build:    17062821
% Date:     Jun-28 2017, 9:40 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2014, 2016, 2017, Jochen Weber
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

% global variable
global ne_gcfg;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% get configuration and handles
cc = ne_gcfg.fcfg;
ch = ne_gcfg.h;

% only works reasonably for VMP files
stvar = cc.StatsVar;
stvix = cc.StatsVarIdx;
if isempty(stvix)
    stvix = 1;
else
    stvix = stvix(1);
end
if ~isxff(stvar, 'vmp')
    return;
end

% remote mode
if nargin > 2 && ischar(varargin{3}) && strcmpi(varargin{3}(:)', 'remote')
    remote = true;
else
    remote = false;
end

% allow only one instance
if any(strcmp('alphasim', ne_gcfg.c.blockcb))
    return;
end
ne_gcfg.c.blockcb{end+1} = 'alphasim';

% get preliminary info
res = stvar.Resolution;
mmap = stvar.Map(stvix);
ntails = 2 - double(mmap.Type == 4);
utails = min(ntails, max(1, mmap.ShowPositiveNegativeFlag - 1));

% compute basic null-distribution shift
zshift = 0;
if any(mmap.Type == [1, 2])
    mstat = double(mmap.VMPData(:));
    mstat(isinf(mstat) | isnan(mstat) | mstat == 0) = [];
    switch mmap.Type

        % t-Map
        case {1}

            % compute z-stat
            z = -sign(mstat) .* sdist('norminv', sdist('tcdf', -abs(mstat), mmap.DF1), 0, 1);

        % r-Map
        case {2}

            % compute z-stat
            z = -sign(mstat) .* sdist('norminv', correlpvalue(abs(mstat), mmap.DF1 + 2), 0, 1);
    end

    % remove voxels we (tend to) believe are not global signal artefacts
    z(abs(z) > 5) = [];

    % detect shift
    zsiter = 0;
    zlshift = 0;
    zshift = median(z(z > (zlshift - 1) & z < (zlshift + 1)));
    while zsiter < 20 && abs(zlshift - zshift) > 0.01
        zsiter = zsiter + 1;
        zlshift = zshift;
        zshift = median(z(z > (zlshift - 1) & z < (zlshift + 1)));
    end
end

% get smoothing kernel size in mm
smkest = ' 12';
writeres = false;
showwarn = true;
if isfield(stvar.Map, 'RunTimeVars') && isstruct(mmap.RunTimeVars) && ...
    isfield(mmap.RunTimeVars, 'FWHMResEst') && ...
    isa(mmap.RunTimeVars.FWHMResEst, 'double') && ...
    numel(mmap.RunTimeVars.FWHMResEst) == 3 && ...
   ~any(isinf(mmap.RunTimeVars.FWHMResEst) | isnan(mmap.RunTimeVars.FWHMResEst))
    smkest = sprintf(' %.1f', harmmean(mmap.RunTimeVars.FWHMResEst));
    writeres = true;
elseif isfield(stvar.Map, 'RunTimeVars') && isstruct(mmap.RunTimeVars) && ...
    isfield(mmap.RunTimeVars, 'FWHMMapEst') && ...
    isa(mmap.RunTimeVars.FWHMMapEst, 'double') && ...
    numel(mmap.RunTimeVars.FWHMMapEst) == 3 && ...
   ~any(isinf(mmap.RunTimeVars.FWHMMapEst) | isnan(mmap.RunTimeVars.FWHMMapEst))
    smkest = sprintf(' %.1f', harmmean(mmap.RunTimeVars.FWHMMapEst));
    writeres = true;
else
    if remote
        vanswer = 'Yes';
    else
        vanswer = questdlg( ...
            'Smoothness of map is not established. Would you like to estimate it now?', ...
            'NeuroElf - user input', 'Yes', 'No', 'Yes');
    end
    if ischar(vanswer) && ~isempty(vanswer) && strcmp(vanswer, 'Yes')
        if ~isfield(stvar.Map, 'RunTimeVars') || ...
           ~isstruct(mmap.RunTimeVars) || numel(mmap.RunTimeVars) ~= 1
            stvar.Map(stvix).RunTimeVars = struct;
        end
        try
            dmap = double(mmap.VMPData);
            dmap(isinf(dmap) | isnan(dmap)) = 0;
            stvar.Map(stvix).RunTimeVars.FWHMMapEst = mapestsmooth(dmap, res);
        catch ne_eo;
            uiwait(errordlg(ne_eo.message, 'NeuroElf - error', 'modal'));
            ne_gcfg.c.blockcb(strcmp('alphasim', ne_gcfg.c.blockcb)) = [];
            return;
        end
        stvar.RunTimeVars.AutoSave = true;
        smkest = sprintf(' %.1f', harmmean(stvar.Map(stvix).RunTimeVars.FWHMMapEst));
        mmap = stvar.Map(stvix);
        writeres = true;
    end
end
if remote
    uinput = {smkest, ' 1000', ' 0.05', ' 0', ' 1'};
    for vc = 4:(min(7, nargin))
        if isa(varargin{vc}, 'double') && numel(varargin{vc}) == 1 && ...
           ~isinf(varargin{vc}) && ~isnan(varargin{vc}) && varargin{vc} >= 0
            uinput{vc-2} = sprintf('%g', varargin{vc});
        end
    end
    showwarn = false;
else
    uinput = inputdlg( ...
       {'Smoothing kernel size (mm):', 'Number of iterations:', 'FWE p-level:', ...
        sprintf('Shift in Z-distribution (detected: %.2f)', zshift), ...
        'Conjunction of how many tests:'}, ...
       'NeuroElf GUI user input', 1, {smkest, ' 1000', ' 0.05', ' 0', ' 1'});
end
if ~iscell(uinput) || numel(uinput) ~= 5 || ~ischar(uinput{1}) || isempty(uinput{1}) || ...
   ~ischar(uinput{2}) || isempty(uinput{2}) || ~ischar(uinput{3}) || isempty(uinput{3}) || ...
   ~ischar(uinput{4}) || isempty(uinput{4}) || ~ischar(uinput{5}) || isempty(uinput{5})
    ne_gcfg.c.blockcb(strcmp('alphasim', ne_gcfg.c.blockcb)) = [];
    return;
end
smkmm = uinput{1};
try
    smkmm = str2double(smkmm);
catch ne_eo;
    warning('neuroelf:general:badInput', ...
        'Invalid user input (kernel size): %s.', ne_eo.message);
    ne_gcfg.c.blockcb(strcmp('alphasim', ne_gcfg.c.blockcb)) = [];
    return;
end

% get number of iterations
niter = uinput{2};
try
    niter = round(str2double(niter));
catch ne_eo;
    warning('neuroelf:general:badInput', ...
        'Invalid user input (number of iterations): %s.', ne_eo.message);
    ne_gcfg.c.blockcb(strcmp('alphasim', ne_gcfg.c.blockcb)) = [];
    return;
end

% get FWE p-level
fwep = uinput{3};
try
    fwep = str2double(fwep);
    if numel(fwep) ~= 1 || isinf(fwep) || isnan(fwep) || fwep <= 0 || fwep > 0.5
        error('Bad number entered.');
    end
catch ne_eo;
    warning('neuroelf:general:badInput', ...
        'Invalid user input (fwe p-level): %s.', ne_eo.message);
end

% get z-shift
zshift = uinput{4};
try
    zshift = str2double(zshift);
    if numel(zshift) ~= 1 || isinf(zshift) || isnan(zshift) || abs(zshift) > 5
        error('Bad number entered.');
    end
catch ne_eo;
    warning('neuroelf:general:badInput', ...
        'Invalid user input (z-shift): %s.', ne_eo.message);
end

% get number of tests in conjunction
nconj = uinput{5};
try
    nconj = str2double(nconj);
    if numel(nconj) ~= 1 || isinf(nconj) || isnan(nconj) || ~any((1:5) == nconj)
        error('Bad number entered.');
    end
catch ne_eo;
    warning('neuroelf:general:badInput', ...
        'Invalid user input (noconj): %s.', ne_eo.message);
end

% get map size and create mask
try

    % take full map (with all non-zero values)
    mvol = (mmap.VMPData ~= 0 & ~isinfnan(mmap.VMPData));

    % check mask
    if sum(mvol(:)) < 16
        mvol = [];
    end

% if an error occurred
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;

    % don't apply mask
    mvol = [];
end

% warn user about issues with alphasim (see Eklund et al. 2016)
% http://www.pnas.org/content/113/28/7900.short
if showwarn
    uiwait(warndlg([ ...
        'Please be advised that, for relatively less stringent thresholds,' char(10) ...
        'such as p < 0.005 or higher, alphasim currently does not adequately' char(10) ...
        'control for the false-positive rate, given the non-Gaussian shape' char(10) ...
        'of the spatial auto-correlation function of the residual, as well' char(10) ...
        'as the non-uniform distribution of smoothness across the brain.' char(10) char(10) ...
        'See http://www.pnas.org/content/113/28/7900.short (Eklund et al. 2016)' char(10) ...
        'for further information!' char(10) char(10) ...
        'As an alternative, please consider using non-parametric cluster-size' char(10) ...
        'correction methods (provided in the contrast manager).'], 'NeuroElf - info', 'modal'));
end

% use the list of thresholds in figure config
alphasthr = ne_gcfg.fcfg.asimthr;

% add additional threshold to be tested
if any(mmap.Type == [1, 2, 4])
    switch (mmap.Type)
        case {1}
            cpthr = sdist('tcdf', -abs(mmap.LowerThreshold), mmap.DF1);
            if utails == 2
                cpthr = 2 * cpthr;
            end
        case {2}
            cpthr = correlpvalue(abs(mmap.LowerThreshold), mmap.DF1 + 2, utails == 1);
        case {4}
            cpthr = sdist('fcdf', abs(mmap.LowerThreshold), mmap.DF1, mmap.DF2, true);
    end
    if ~any(alphasthr > (0.999 * cpthr) & alphasthr < (1.001 * cpthr))
        alphasthr(end+1) = cpthr;
    end
end

% make progress bar visible
cprog = ne_progress(0, 0, {true, 0, 'alphasim'});

% echo
if ne_gcfg.c.echo
    ne_echo({['alphasim_results = alphasim(%s, struct(' ...
        '''conj'', %d, ''mask'', maskvol, ''fwhm'', %s, ''niter'', %d, ' ...
        '''stype'', [%d, %d], ''thr'', %s, ''zshift'', %s);'], ...
        nconj, any2ascii(size(mvol)), any2ascii(smkmm([1, 1, 1]) ./ res), ...
        niter, utails, ntails, any2ascii(alphasthr), any2ascii(zshift)});
end

% run alphasim
alphasout = alphasim(size(mvol), struct('conj', nconj, 'mask', mvol, ...
    'fwhm', (1 / res) .* smkmm([1, 1, 1]), 'niter', niter, ...
    'stype', [utails, ntails], 'thr', alphasthr, ...
    'zshift', zshift, 'pbar', ch.Progress));

% then hide progress bar again
ne_progress(0, 0, cprog);

% output result in console
disp(' ');
disp('alphasim output for');
fprintf(' - %d-by-%d-by-%d volume with %d voxels of %dmm resolution\n', ...
    [size(mvol), sum(mvol(:)), res]);
fprintf(' - %.1fmm smoothing kernel\n', smkmm);
fprintf(' - %d iterations\n', niter);
fprintf(' - %d of %d tails\n', utails, ntails);
disp(repmat('-', 1, 72));

% for each threshold
alphakthr = alphasthr;
for tc = 1:numel(alphasthr)

    % find first occurrence with output p < FWE-level
    kthr = findfirst(alphasout{tc}(:,end) < fwep);

    % if isempty, last size didn't succeed yet
    if isempty(kthr)
        kthr = size(alphasout{tc}, 1) + 1;
    end

    % store in list
    alphakthr(tc) = kthr;

    % and display line
    fprintf('FWE(p < %g) k for uncorr. p < %g : %d voxels\n', ...
        fwep, alphasthr(tc), kthr);
end
disp(' ');

% if current object is VMP and single map and clustercheck is enabled
if numel(cc.StatsVarIdx) == 1 && ch.Stats.UsekThr.Value > 0

    % figure out if we can usefully update cluster threshold
    try

        % get threshold from dropdown
        setthr = str2double(ch.Stats.PThresh.String{ch.Stats.PThresh.Value});

        % if any value is a match
        thrmatch = (alphasthr == setthr);
        if any(thrmatch)

            % then update cluster size
            set(ch.Stats.kThresh, 'String', sprintf('%d', alphakthr(thrmatch)));

            % and run the corresponding callback
            ne_setstatthrccheck;
        end
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
end

% store results
if writeres
    fld = sprintf('AlphaSim%d%d', utails, ntails);
    stvar.Map(stvix).RunTimeVars.(fld) = {fwep, [alphasthr(:), alphakthr(:)]};
end

% then unblock another run
ne_gcfg.c.blockcb(strcmp('alphasim', ne_gcfg.c.blockcb)) = [];
