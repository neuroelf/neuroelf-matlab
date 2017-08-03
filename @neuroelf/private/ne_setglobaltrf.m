% FUNCTION ne_setglobaltrf: set global transformation
function ne_setglobaltrf(varargin)

% Version:  v0.9c
% Build:    11052400
% Date:     Apr-29 2011, 8:11 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, Jochen Weber
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

% position of controls on which clicks allows updating of position
cc = ne_gcfg.fcfg;
ch = ne_gcfg.h;
fPos = [cc.slicepos; cc.zslicepos];
nPos = ch.MainFig.CurrentPoint;

% do we need a hit-test
if any(ne_gcfg.c.btdown == gcbf) && ...
    isempty(cc.mpos.ddat)

    % for first page (3-slice view)
    if cc.page == 1

        % set object positions of SAG/COR/TRA image to -1 (cannot hit)
        fPos(4:6, :) = -1;

        % make hit-test
        cobj = findfirst( ...
            fPos(:, 1) <= nPos(1) & fPos(:, 2) <= nPos(2) & ...
            fPos(:, 3) >  nPos(1) & fPos(:, 4) >  nPos(2));

    % for page 2
    elseif cc.page == 2

        % do the reverse (set 3-slice objects' positions to -1)
        fPos(1:3, :) = -1;
        cobj = findfirst( ...
            fPos(:, 1) <= nPos(1) & fPos(:, 2) <= nPos(2) & ...
            fPos(:, 3) >  nPos(1) & fPos(:, 4) >  nPos(2));
        if ~isempty(cobj)
            cobj = cobj + cc.zoom - 1;
        end
    else
        cobj = [];
    end

    % object not hit?
    if isempty(cobj)

        % then pretend the button wasn't hit after all
        ne_gcfg.c.btdown = [];
        return;
    end

    % get config from current variable
    slvar = cc.SliceVar;
    if numel(slvar) == 1 && ...
        isxff(slvar, true) && ...
        isfield(slvar.RunTimeVars, 'TrfPlus') && ...
        isequal(size(slvar.RunTimeVars.TrfPlus), [4, 4])
        sltrf = spmitrf(slvar.RunTimeVars.TrfPlus);
    else
        sltrf = spmitrf(eye(4));
    end

    % store config
    ne_gcfg.fcfg.mpos.ddat = {cc.page, nPos, cc, cobj, sltrf};

% we still need to update the position
elseif any(ne_gcfg.c.btdown == gcbf)

    % get original position and configuration
    mods = cc.mpos.mods;
    sltrf = cc.mpos.ddat{5};
    oPos = cc.mpos.ddat{2};
    cobj = cc.mpos.ddat{4};
    cc = cc.mpos.ddat{3};

    % shift -> translation
    if numel(mods) == 1 && ...
        strcmpi(mods{1}, 'shift')

        % compute new translation
        switch (cobj)
            case {1}
                cc.strtra(2) = min(128, max(-128, ...
                    cc.strtra(2) + nPos(1) - oPos(1)));
                cc.strtra(3) = min(128, max(-128, ...
                    cc.strtra(3) + oPos(2) - nPos(2)));
            case {2}
                cc.strtra(1) = min(128, max(-128, ...
                    cc.strtra(1) + oPos(1) - nPos(1)));
                cc.strtra(3) = min(128, max(-128, ...
                    cc.strtra(3) + oPos(2) - nPos(2)));
            case {3}
                cc.strtra(1) = min(128, max(-128, ...
                    cc.strtra(1) + oPos(1) - nPos(1)));
                cc.strtra(2) = min(128, max(-128, ...
                    cc.strtra(2) + oPos(2) - nPos(2)));
            case {4}
                cc.strtra(2) = min(128, max(-128, ...
                    cc.strtra(2) + 0.5 * (nPos(1) - oPos(1))));
                cc.strtra(3) = min(128, max(-128, ...
                    cc.strtra(3) + 0.5 * (oPos(2) - nPos(2))));
            case {5}
                cc.strtra(1) = min(128, max(-128, ...
                    cc.strtra(1) + 0.5 * (oPos(1) - nPos(1))));
                cc.strtra(3) = min(128, max(-128, ...
                    cc.strtra(3) + 0.5 * (oPos(2) - nPos(2))));
            case {6}
                cc.strtra(1) = min(128, max(-128, ...
                    cc.strtra(1) + 0.5 * (oPos(1) - nPos(1))));
                cc.strtra(2) = min(128, max(-128, ...
                    cc.strtra(2) + 0.5 * (oPos(2) - nPos(2))));
        end

    % alt -> zooming
    elseif numel(mods) == 1 && ...
        strcmpi(mods{1}, 'alt')

        % compute new translation
        switch (cobj)
            case {1, 2, 3}
                cc.strscl = cc.strscl .* (1.01 ^ (round(oPos(2) - nPos(2))));
            case {4, 5, 6}
                cc.strscl = cc.strscl .* (1.005 ^ (round(oPos(2) - nPos(2))));
        end
        if min(abs(cc.strscl)) < 0.1
            cc.strscl = cc.strscl .* (0.1 / min(abs(cc.strscl)));
        end
        if max(abs(cc.strscl)) > 8
            cc.strscl = cc.strscl .* (8 / max(abs(cc.strscl)));
        end

    % alt + shift -> rotation
    elseif numel(mods) == 2 && ...
        any(strcmpi(mods, 'alt')) && ...
        any(strcmpi(mods, 'shift'))

        % compute new rotation
        switch (cobj)
            case {1, 4}
                cc.strrot(1) = mod(cc.strrot(1) + (nPos(1) - oPos(1)), 360);
            case {2, 5}
                cc.strrot(2) = mod(cc.strrot(2) + (oPos(1) - nPos(1)), 360);
            case {3, 6}
                cc.strrot(3) = mod(cc.strrot(3) + (oPos(1) - nPos(1)), 360);
        end

    % crtl + shift -> translation of object
    elseif numel(mods) == 2 && ...
        any(strcmpi(mods, 'control')) && ...
        any(strcmpi(mods, 'shift'))

        % compute new translation
        switch (cobj)
            case {1}
                sltrf{1}(2) = min(256, max(-256, ...
                    sltrf{1}(2) + oPos(1) - nPos(1)));
                sltrf{1}(3) = min(256, max(-256, ...
                    sltrf{1}(3) + nPos(2) - oPos(2)));
            case {2}
                sltrf{1}(1) = min(256, max(-256, ...
                    sltrf{1}(1) + nPos(1) - oPos(1)));
                sltrf{1}(3) = min(256, max(-256, ...
                    sltrf{1}(3) + nPos(2) - oPos(2)));
            case {3}
                sltrf{1}(1) = min(256, max(-256, ...
                    sltrf{1}(1) + nPos(1) - oPos(1)));
                sltrf{1}(2) = min(256, max(-256, ...
                    sltrf{1}(2) + nPos(2) - oPos(2)));
            case {4}
                sltrf{1}(2) = min(256, max(-256, ...
                    sltrf{1}(2) + 0.5 * (oPos(1) - nPos(1))));
                sltrf{1}(3) = min(256, max(-256, ...
                    sltrf{1}(3) + 0.5 * (nPos(2) - oPos(2))));
            case {5}
                sltrf{1}(1) = min(256, max(-256, ...
                    sltrf{1}(1) + 0.5 * (nPos(1) - oPos(1))));
                sltrf{1}(3) = min(256, max(-256, ...
                    sltrf{1}(3) + 0.5 * (nPos(2) - oPos(2))));
            case {6}
                sltrf{1}(1) = min(256, max(-256, ...
                    sltrf{1}(1) + 0.5 * (nPos(1) - oPos(1))));
                sltrf{1}(2) = min(256, max(-256, ...
                    sltrf{1}(2) + 0.5 * (nPos(2) - oPos(2))));
        end

    % crtl + alt -> zooming of object
    elseif numel(mods) == 2 && ...
        any(strcmpi(mods, 'control')) && ...
        any(strcmpi(mods, 'alt'))

        % compute new translation
        switch (cobj)
            case {1, 2, 3}
                sltrf{3} = sltrf{3} .* (1.01 ^ (round(oPos(2) - nPos(2))));
            case {4, 5, 6}
                sltrf{3} = sltrf{3} .* (1.005 ^ (round(oPos(2) - nPos(2))));
        end
        if min(abs(sltrf{3})) < 0.1
            sltrf{3} = sltrf{3} .* (0.1 / min(abs(sltrf{3})));
        end
        if max(abs(sltrf{3})) > 8
            sltrf{3} = sltrf{3} .* (8 / max(abs(sltrf{3})));
        end

    % crtr + alt + shift -> rotation
    elseif numel(mods) == 3 && ...
        any(strcmpi(mods, 'control')) && ...
        any(strcmpi(mods, 'alt')) && ...
        any(strcmpi(mods, 'shift'))

        % compute new rotation
        pfac = pi / 180;
        switch (cobj)
            case {1, 4}
                sltrf{2}(1) = sltrf{2}(1) + pfac * (oPos(1) - nPos(1));
            case {2, 5}
                sltrf{2}(2) = sltrf{2}(2) + pfac * (nPos(1) - oPos(1));
            case {3, 6}
                sltrf{2}(3) = sltrf{2}(3) + pfac * (nPos(1) - oPos(1));
        end
        sltrfr = sltrf{2} < 0;
        if any(sltrfr)
            sltrf{2}(sltrfr) = sltrf{2}(sltrfr) + ...
                (2 * pi) .* ceil(abs(sltrf{2}(sltrfr)) ./ (2 * pi));
        end
        sltrfr = sltrf{2} >= (2 * pi);
        if any(sltrfr)
            sltrf{2}(sltrfr) = sltrf{2}(sltrfr) - ...
                (2 * pi) .* floor(sltrf{2}(sltrfr) ./ (2 * pi));
        end
    end

    % re-store config
    ne_gcfg.fcfg.strrot = cc.strrot;
    ne_gcfg.fcfg.strtra = cc.strtra;
    ne_gcfg.fcfg.strscl = cc.strscl;
    ne_gcfg.fcfg.strans = spmtrf(cc.strtra, (pi / 180) .* cc.strrot, 1 ./ cc.strscl);
    ne_gcfg.fcfg.SliceVar.RunTimeVars.TrfPlus = spmtrf(sltrf{:});

    % set position
    ne_setslicepos(0, 0, cc.cpos);
end
