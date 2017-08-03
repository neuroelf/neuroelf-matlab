function xo = vmp_Update(xo, F, S, V)
% VMP::Update  - called after subsasgn for VMPs
%
% Using: applyfdr, flexinterpn, makelabel.

% Version:  v1.1
% Build:    16021413
% Date:     Feb-14 2016, 1:03 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2011, 2014, 2016, Jochen Weber
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

% global settings
global xffsngl;

% persistent config
persistent xffvmp;
if isempty(xffvmp)

    % from neuroelf library
    using(neuroelf, {'applyfdr', 'flexinterpn', 'makelabel'});

    % initialize struct
    xffvmp = struct;

    % settings
    xffvmp.applyfdr = applyfdr;
    xffvmp.fdr = xffsngl.CONF.settings.Statistics.FDR.Thresholds;
    xffvmp.flexinterpn = flexinterpn;
    xffvmp.makelabel = makelabel;
end

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'vmp') || ~ischar(F) || isempty(F)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
if nargin < 3 || ~isstruct(S)
    S = struct;
    S.type = '.';
    S.subs = F;
end
if nargin < 4
    V = [];
end

% get content
bc = xo.C;

% linearize
F = xffvmp.makelabel(F(:)');

% F valid
if ~isfield(bc, F)
    error('neuroelf:xff:invalidProperty', 'Cannot find property ''%s'' for type VMP.', F);
end

% only relevant for fileversion
if ~strcmpi(F, 'fileversion')
    return;
end

% get wanted and old version
reqv = bc.FileVersion;
oldv = V;

% check version
if isempty(reqv) || ~isa(reqv, 'double') || isnan(reqv(1)) || isinf(reqv(1)) || ...
    fix(reqv(1)) ~= reqv(1) || reqv(1) < 2 || reqv(1) > 6 || isempty(oldv) || ...
   ~isa(oldv, 'double') || isnan(oldv(1)) || isinf(oldv(1)) || ...
    fix(oldv(1)) ~= oldv(1) || oldv(1) < 2 || oldv(1) > 6

    % give warning if required
    warning('neuroelf:xff:invalidPropertyValue', 'Invalid old/new FileVersion value.');

    % guess better values
    reqv = 5;
    bc.FileVersion = reqv;
    if bc.NativeResolutionFile
        oldv = 5;
    else
        oldv = 4;
    end
end

% if resolution is 1 make sure to go to old format
if bc.Resolution == 1
    reqv = 4;
end

% something to do?
if (oldv < 5 && reqv < 5) || (oldv > 4 && reqv > 4)
    xo.C = bc;
    return;
end

% build new content
newCONT = xffnewcont('vmp');
oldMap  = bc.Map;
newMap  = newCONT.Map;

% upgrade
if oldv < 5

    % set global values
    newCONT.NativeResolutionFile = 1;
    newCONT.NrOfMaps   = numel(oldMap);
    newCONT.XStart     = bc.XStart;
    newCONT.XEnd       = bc.XEnd;
    newCONT.YStart     = bc.YStart;
    newCONT.YEnd       = bc.YEnd;
    newCONT.ZStart     = bc.ZStart;
    newCONT.ZEnd       = bc.ZEnd;
    newCONT.Resolution = bc.Resolution;
    dimres = newCONT.Resolution;
    dimX = (newCONT.XEnd - newCONT.XStart) / dimres;
    dimY = (newCONT.YEnd - newCONT.YStart) / dimres;
    dimZ = (newCONT.ZEnd - newCONT.ZStart) / dimres;
    if any([dimX, dimY, dimZ] ~= fix([dimX, dimY, dimZ]))
        error('neuroelf:xff:dimensionError', 'X/Y/Z Start/End values must match with Resolution.');
    end

    if newCONT.NrOfMaps > 0
        newres = newCONT.Resolution;
        oldres = [(newCONT.XEnd - newCONT.XStart) / size(oldMap(1).VMPData, 1), ...
            (newCONT.YEnd - newCONT.YStart) / size(oldMap(1).VMPData, 2), ...
            (newCONT.ZEnd - newCONT.ZStart) / size(oldMap(1).VMPData, 3)];
        if any(oldres ~= fix(oldres))
            oldres = [(1 + newCONT.XEnd - newCONT.XStart) / size(oldMap(1).VMPData, 1), ...
                (1 + newCONT.YEnd - newCONT.YStart) / size(oldMap(1).VMPData, 2), ...
                (1 + newCONT.ZEnd - newCONT.ZStart) / size(oldMap(1).VMPData, 3)];
        end
        if any(diff(oldres)) || any(oldres ~= fix(oldres))
            error('neuroelf:xff:dimensionError', 'Map dimension mismatch with in X/Y/Z Start/End.');
        end
        oldres = oldres(1);
        difres = newres / oldres;
        if difres ~= fix(difres)
            scrds = 0.01 + (1:difres:(difres * (max([dimX, dimY, dimZ]) + 1)));
        end
    end

    % iterate over maps
    for mc = 1:numel(oldMap)

        % set values
        newMap(mc).Type               = oldMap(mc).Type;
        newMap(mc).LowerThreshold     = oldMap(mc).LowerThreshold;
        newMap(mc).UpperThreshold     = oldMap(mc).UpperThreshold;
        newMap(mc).Name               = oldMap(mc).Name;
        newMap(mc).RGBLowerThreshPos  = oldMap(mc).RGBLowerThreshPos;
        newMap(mc).RGBUpperThreshPos  = oldMap(mc).RGBUpperThreshPos;
        newMap(mc).RGBLowerThreshNeg  = oldMap(mc).RGBLowerThreshNeg;
        newMap(mc).RGBUpperThreshNeg  = oldMap(mc).RGBUpperThreshNeg;
        newMap(mc).UseRGBColor        = oldMap(mc).UseRGBColor;
        newMap(mc).LUTName            = oldMap(mc).LUTName;
        newMap(mc).TransColorFactor   = oldMap(mc).TransColorFactor;
        newMap(mc).NrOfLags           = oldMap(mc).NrOfLags;
        newMap(mc).MinLag             = oldMap(mc).MinLag;
        newMap(mc).MaxLag             = oldMap(mc).MaxLag;
        newMap(mc).CCOverlay          = oldMap(mc).CCOverlay;
        newMap(mc).ClusterSize        = ...
            round(oldMap(mc).ClusterSize / difres ^ 3);
        newMap(mc).EnableClusterCheck = oldMap(mc).EnableClusterCheck;
        newMap(mc).UseValuesAboveThresh = ...
            oldMap(mc).UseValuesAboveThresh;
        newMap(mc).DF1                = oldMap(mc).DF1;
        newMap(mc).DF2                = oldMap(mc).DF2;
        newMap(mc).ShowPositiveNegativeFlag = 3;
        newMap(mc).BonferroniValue    = oldMap(mc).BonferroniValue;
        newMap(mc).NrOfFDRThresholds  = numel(xffvmp.fdr);
        newMap(mc).UnknownValue       = -1;
        newMap(mc).TimePointData      = zeros(0, 1);
        if difres == fix(difres)
            newMap(mc).VMPData        = ...
                oldMap(mc).VMPData(1:difres:end-1, 1:difres:end-1, 1:difres:end-1);
        else
            newMap(mc).VMPData        = ...
                xffvmp.flexinterpn(oldMap(mc).VMPData, ...
                    [Inf, Inf, Inf; 1, 1, 1; difres, difres, difres; ...
                     scrds(dimX), scrds(dimY), scrds(dimZ)], [0;1;0], 1, 0);
        end

        % build FDR table
        try
            mvals = newMap(mc).VMPData(~isnan(newMap(mc).VMPData(:)) & ...
                (newMap(mc).VMPData(:) ~= 0));
            if ~isempty(mvals)
                switch (newMap(mc).Type)
                    case 1  % t-score
                        newMap(mc).FDRThresholds = ...
                            [xffvmp.fdr(:), ...
                            xffvmp.applyfdr(double(mvals), 't', ...
                            xffvmp.fdr(:), newMap(mc).DF1, ...
                            [], true)];
                    case 2  % correlation
                        newMap(mc).FDRThresholds = ...
                            [xffvmp.fdr(:), ...
                            xffvmp.applyfdr(double(mvals), 'r', ...
                            xffvmp.fdr(:), newMap(mc).DF1, ...
                            [], true)];
                    case 4  % F-score
                        newMap(mc).FDRThresholds = ...
                            [xffvmp.fdr(:), ...
                            xffvmp.applyfdr(double(mvals), 'F', ...
                            xffvmp.fdr(:), newMap(mc).DF1, ...
                            newMap(mc).DF2, true)];
                    otherwise
                        error('FDR_ERROR');
                end
            else
                newMap(mc).NrOfFDRThresholds = 1;
                newMap(mc).FDRThresholds = [0, 1e5, 1e5];
            end
        catch xfferror
            neuroelf_lasterr(xfferror);
            newMap(mc).NrOfFDRThresholds = 1;
            newMap(mc).FDRThresholds = [0, 1e5, 1e5];
        end
    end

    % put Map into new object
    newCONT.Map = newMap;

% downgrade
else
    try
        nxo = vmp_MakeHiResVMP(xo);
    catch xfferror
        warning('neuroelf:xff:conversionError', ...
            'Error converting NatRes to HiRes VMP: %s.', xfferror.message);
        return;
    end
    newCONT = nxo.C;
    delete(nxo);
    newCONT.FileVersion = reqv;
end

% put back into object
bc = newCONT;

% set back
xo.C = bc;
