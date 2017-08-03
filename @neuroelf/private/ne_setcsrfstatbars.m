function varargout = ne_setcsrfstatbars(varargin)
% ne_setcsrfstatbars  - set current Surfaces stats color bars
%
% FORMAT:       ne_setcsrfstatbars(SRC, EVT [, window [, showbars]])
%
% Input fields:
%
%       SRC, EVT    Matlab handle callback inputs (discarded)
%       window      window specifier (used to target satellite windows)
%       showbars    override global flag (either 'on' or 'off')
%
% No output fields.
%
% Example:
%
%       [surface_viewer, ~, surface_id] = neuroelf_gui('undock');
%       neuroelf_gui('satresize', surface_id, [640, 480]);
%       ne_setcsrfstatbars(0, 0, surface_id, 'off');

% Version:  v1.1
% Build:    16052922
% Date:     May-29 2016, 10:02 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2014, 2016, Jochen Weber
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

% which window
if nargin < 3 || ~ischar(varargin{3}) || isempty(varargin{3}) || ...
   ~isfield(ne_gcfg.cc, varargin{3}(:)')
    cc = ne_gcfg.fcfg;
    ch = ne_gcfg.h;
    scn = mlhandle(ch.Scenery);
    sci = get(scn, 'Value');
    scu = get(scn, 'UserData');
else
    ch = ne_gcfg.cc.(varargin{3}(:)');
    cc = ch.Config;
    scn = ch.Scenery;
    sci = scn.Value;
    scu = scn.UserData;
end
cinist = ne_gcfg.c.ini.Statistics;
if nargin < 4 || ~ischar(varargin{4}) || isempty(varargin{4}) || ...
   ~any(strcmpi(varargin{4}(:)', {'on', 'off'}))
    showbars = cinist.ShowThreshBars;
else
    showbars = strcmpi(varargin{4}(:)', 'on');
end
showtext = cinist.ShowThreshText;
sdtxt = sprintf('%%.%df', cinist.ShowThreshDecimals);

% don't show stats bars?
set(ch.SurfaceStatsText, 'Visible', 'off');
if ~showbars || isempty(scu) || isempty(sci) || cc.page ~= 3
    set(ch.SurfaceStatsBar, 'Visible', 'off');
    return;
end

% start with scenery
sbars = cell(1, numel(sci));
sbarl = sbars;
sbaru = sbars;
sbart = sbars;
for uc = 1:numel(sci)
    f = scu{sci(uc), 4};
    try
        fh = handles(f);
        sbars{uc} = fh.SurfStatsBars;
        if isxff(fh.Stats{1}, {'fsmf', 'mtc', 'smp'}) && ~isempty(fh.Stats{2})
            sbarms = fh.Stats{1}.Map(fh.Stats{2});
            sbarul = zeros(2, numel(fh.Stats{2}));
            sbarut = zeros(3, size(sbarul, 2));
            for umc = 1:size(sbarul, 2)
                sbarul(1, umc) = sbarms(umc).LowerThreshold;
                sbarul(2, umc) = sbarms(umc).UpperThreshold;
                sbarut(:, umc) = [sbarms(umc).Type; sbarms(umc).DF1; sbarms(umc).DF2];
            end
            sbarl{uc} = sbarul(1, :);
            sbaru{uc} = sbarul(2, :);
            sbart{uc} = sbarut;
        end
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
end
sbarsz = cc.SurfBarSize;
sbars = cat(2, sbars{:});
sbarv = cat(2, sbars{:});
sbarw = size(sbarv, 2);
if sbarw == 0
    vflag = get(ch.SurfaceStatsBar, 'Visible');
    if ~strcmpi(vflag, 'off')
        set(ch.SurfaceStatsBar, 'Visible', 'off');
    end
    return;
end
sbarf = (sbarsz(2) +1 - sbarw) / sbarw;
blanks = 1;
if sbarf < 1
    blanks = 0;
    sbarf = floor(sbarsz(2) / sbarw);
    if sbarf < 1
        return;
    end
end
sbari = floor(1:(sbarf+blanks):(sbarsz(2)+0.9));
sbari(end+1) = sbarsz(2) + 1;
sbari = sbari(:)';

% create image first
barimage = repmat(reshape(cc.SurfBackColor, [1, 1, 3]), sbarsz);
for bc = 1:sbarw
    sbarf = sbari(bc+1)-(sbari(bc) + blanks);
    if sbarf > 1
        barimage(:, sbari(bc):sbari(bc)+sbarf-1, :) = repmat(sbarv(:, bc, :), 1, sbarf);
    else
        barimage(:, sbari(bc), :) = sbarv(:, bc, :);
    end
end

% output image?
if nargout > 0
    varargout{1} = barimage;
end

% prepare for display
barimage = reshape(barimage, [1, prod(sbarsz), 3]);
barimage = reshape(cat(1, barimage, barimage), [2 * prod(sbarsz), 3]);

% set bar visible if data
if sbarw > 0
    set(ch.SurfaceStatsBar, 'FaceVertexCData', barimage, 'Visible', 'on');

    % and check the stats
    sbarl = cat(2, sbarl{:});
    sbaru = cat(2, sbaru{:});
    sbart = cat(2, sbart{:});
    if all(sbarl == sbarl(1)) && all(sbaru == sbaru(1)) && ...
        isequal(sbart, repmat(sbart(:, 1), 1, size(sbart, 2)))
        sbardf1 = sbart(2, 1);
        sbardf2 = sbart(3, 1);
        sbart = sbart(1);
        switch (sbart)
            case 1
                sbartu = sprintf(['t_{[%d]}<' sdtxt], sbardf1, sbaru(1));
                sbart = sprintf(['t_{[%d]}>=' sdtxt], sbardf1, sbarl(1));
            case 2
                sbartu = sprintf(['r_{df=%d}<' sdtxt], sbardf1, sbaru(1));
                sbart = sprintf(['r_{df=%d}>=' sdtxt], sbardf1, sbarl(1));
            case 4
                sbartu = sprintf(['F_{[%d,%d]}<' sdtxt], sbardf1, sbardf2, sbaru(1));
                sbart = sprintf(['F_{[%d,%d]}>=' sdtxt], sbardf1, sbardf2, sbarl(1));
            otherwise
                sbart = '';
        end
        if ~isempty(sbart)
            set(ch.SurfaceStatsText(1), 'String', sbart);
            set(ch.SurfaceStatsText(2), 'String', sbartu);
            if showtext
                set(ch.SurfaceStatsText, 'Visible', 'on');
            end
        end
    end
else
    set(ch.SurfaceStatsBar, 'Visible', 'off');
end
