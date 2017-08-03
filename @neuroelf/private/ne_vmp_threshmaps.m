% FUNCTION ne_vmp_threshmaps: threshold several maps
function ne_vmp_threshmaps(varargin)

% Version:  v1.1
% Build:    16031016
% Date:     Mar-10 2016, 4:34 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011 - 2016, Jochen Weber
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

% only valid if single VMP
stvar = ne_gcfg.fcfg.StatsVar;
stv = true;
if numel(stvar) ~= 1 || ...
   ~isxff(stvar, {'glm', 'vmp'})
    stv = false;
    stvar = ne_gcfg.fcfg.SurfStatsVar;
    if ~isxff(stvar, {'glm', 'smp'})
        return;
    end
end
if stv
    stvix = ne_gcfg.fcfg.StatsVarIdx;
else
    stvix = ne_gcfg.fcfg.SurfStatsVarIdx;
end
stype = stvar.Filetype;

% show a selector
stvarn = stvar.MapNames(true);
if isempty(stvarn)
    return;
end
[csel, csok] = listdlg( ...
    'PromptString', 'Please select maps to threshold...', ...
    'ListString',   stvarn, ...
    'InitialValue', stvix(:), ...
    'ListSize',     [min(600, max(300, 10 .* size(char(stvarn), 2))), 200]);
if ~isequal(csok, 1) || ...
    isempty(csel)
    return;
end
if isempty(stvix)
    stvix = csel;
end

% get maps
maps = stvar.Map;

% request threshold to use (from first selected map)
if nargin > 2 && ...
    iscell(varargin{3}) && ...
    numel(varargin{3}) && ...
   ~any(cellfun('isempty', varargin{3})) && ...
   ~any(isnan(str2double(varargin{3})))
    thr = varargin{3}(:)';
elseif ~strcmpi(stype, 'glm')
    thr = inputdlg({'Lower threshold value (<= 0.1: probability):', ...
        'Upper threshold value (<= 0.1: probability):', ...
        'Transparency value (2 = leave as it is):'}, ...
        'NeuroElf - user input', 1, ...
        {sprintf('%g', maps(stvix(1)).LowerThreshold), ...
         sprintf('%g', maps(stvix(1)).UpperThreshold), ...
         sprintf('%g', maps(stvix(1)).TransColorFactor)});
else
    thr = inputdlg({'Lower threshold value:', ...
        'Upper threshold value:', 'Transparency value (2 = leave as is)'}, ...
        'NeuroElf - user input', 1, ...
        {sprintf('%g', maps(stvix(1)).LowerThreshold), ...
         sprintf('%g', maps(stvix(1)).UpperThreshold), ...
         sprintf('%g', maps(stvix(1)).TransColorFactor)});
end

% test input
if ~iscell(thr) || ...
    numel(thr) ~= 3 || ...
    isempty(thr{1}) || ...
    isempty(thr{2}) || ...
    isempty(thr{3})
    return;
end
try
    lthr = str2double(thr{1});
    uthr = str2double(thr{2});
    newa = str2double(thr{3});
    if numel(lthr) ~= 1 || ...
        isinf(lthr) || ...
        isnan(lthr) || ...
        lthr <= 0 || ...
        numel(uthr) ~= 1 || ...
        isinf(uthr) || ...
        isnan(uthr) || ...
        uthr <= 0
        return;
    end
    if isinf(newa) || ...
        isnan(newa) || ...
        newa > 1
        newa = 2;
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    return;
end

% see if they are compatible
maps = maps(csel);
mapt = cat(1, maps.Type);
mapd1 = cat(1, maps.DF1);
mapd2 = cat(1, maps.DF2);

% only if DF1 and, if required, DF2 are matching
if all(mapt == mapt(1)) && ...
    all(mapd1 == mapd1(1)) && ...
   (mapt(1) ~= 4 || ...
    all(mapd2 == mapd2(1)))

    % get common thresholds
    [l1thr, l2thr, u1thr, u2thr] = ...
        mapthr(mapt(1), mapd1(1), mapd2(1), lthr, uthr);

    % then assign to maps
    for mc = 1:numel(maps)
        if maps(mc).ShowPositiveNegativeFlag > 2
            maps(mc).LowerThreshold = l2thr;
            maps(mc).UpperThreshold = u2thr;
        else
            maps(mc).LowerThreshold = l1thr;
            maps(mc).UpperThreshold = u1thr;
        end
        if newa ~= 2
            maps(mc).TransColorFactor = newa;
        end
    end

% individual thresholding
else

    % iterate first
    for mc = 1:numel(maps)

        % get thresholds
        [l1thr, l2thr, u1thr, u2thr] = ...
            mapthr(mapt(mc), mapd1(mc), mapd2(mc), lthr, uthr);

        % then assign
        if maps(mc).ShowPositiveNegativeFlag > 2
            maps(mc).LowerThreshold = l2thr;
            maps(mc).UpperThreshold = u2thr;
        else
            maps(mc).LowerThreshold = l1thr;
            maps(mc).UpperThreshold = u1thr;
        end
        if newa ~= 2
            maps(mc).TransColorFactor = newa;
        end
    end
end

% set back
stvar.Map(csel) = maps;

% update
if stv
    ne_openfile(0, 0, stvar);
else
    ne_setcsrfstatmap;
end

% sub-function
function [l1thr, l2thr, u1thr, u2thr] = mapthr(mapt, d1, d2, lthr, uthr)

% default
l1thr = lthr;
l2thr = lthr;
u1thr = uthr;
u2thr = uthr;

% depending on type
switch mapt

    % t-map
    case {1}
        if lthr <= 0.1
            l1thr = -sdist('tinv', lthr, d1);
            l2thr = -sdist('tinv', 0.5 * lthr, d1);
        end
        if uthr <= 0.1
            u1thr = -sdist('tinv', uthr, d1);
            u2thr = -sdist('tinv', 0.5 * uthr, d1);
        end
        if u1thr <= l1thr
            u1thr = l1thr - sdist('tinv', 0.05, d1);
        end
        if u2thr <= l2thr
            u2thr = l2thr - sdist('tinv', 0.025, d1);
        end

    % r-map
    case {2}
        if lthr <= 0.1
            l1thr = correlinvtstat(-sdist('tinv', 2 * lthr, d1), ...
                d1 + 2);
            l2thr = correlinvtstat(-sdist('tinv', lthr, d1), ...
                d1 + 2);
        end
        if uthr <= 0.1
            u1thr = correlinvtstat(-sdist('tinv', 2 * uthr, d1), ...
                d1 + 2);
            u2thr = correlinvtstat(-sdist('tinv', uthr, d1), ...
                d1 + 2);
        end
        if u1thr <= l1thr
            u1thr = 1 - 0.5 * l1thr;
        end
        if u2thr <= l2thr
            u2thr = 1 - 0.5 * l2thr;
        end

    % F-map
    case {4}
        if lthr <= 0.1
            l1thr = sdist('finv', lthr, d1, d2, true);
        end
        if uthr <= 0.1
            u1thr = sdist('finv', uthr, d1, d2, true);
        end
        if u1thr <= l1thr
            u1thr = l1thr * 1.44;
        end
        l2thr = l1thr;
        u2thr = u1thr;

    % z-map
    case {12}
        if lthr <= 0.1
            l1thr = -sdist('tinv', lthr, d1);
            l2thr = -sdist('tinv', 0.5 * lthr, d1);
        end
        if uthr <= 0.1
            u1thr = -sdist('tinv', uthr, d1);
            u2thr = -sdist('tinv', 0.5 * uthr, d1);
        end
        if u1thr <= l1thr
            u1thr = l1thr - sdist('tinv', 0.05, d1);
        end
        if u2thr <= l2thr
            u2thr = l2thr - sdist('tinv', 0.025, d1);
        end

    % otherwise
    otherwise

        % simply use values as given
        if uthr <= lthr
            uthr = 1.25 * lthr;
        end
        u1thr = uthr;
        u2thr = uthr;
end
