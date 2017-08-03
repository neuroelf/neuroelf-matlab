function xo = aft_SetColors(xo, maps, cspec)
% AFT::SetColors  - set the OverlayColors for the GUI
%
% FORMAT:       [obj = ] obj.SetColors([maps, [cspec]])
%
% Input fields:
%
%       maps        optional specification of map number(s)
%       cspec       optional color specification, either of
%                   - Cx3 RGB colors (even number, range [0...255])
%                   - OLT object (from which the .Colors are taken)
%                   - OLT filename (will be loaded temporarily)
%                   - 'auto' apply RGB and LUTs as specified
%                   - 'RGB' to temporarily override Map setting
%                   - 'HSV' convert RGB to HSV and then grade
%                   - 'HSV4RGB' only apply HSV of RGB is selected
%                   - 'xauto', same as auto but use extended for <default>
%                   - 'xrgb', use 40 instead of 10 RGB colors (per tail)
%                   - 'xhsv', use 40 instead of 10 HSV grades (per tail)
%                   - empty array (set OverlayColors to [])
%                   with the default being to 'auto'
%
% TYPES: AVA, CMP, FSMF, GLM, HDR, HEAD, MTC, SMP, VMP, VTC

% Version:  v1.1
% Build:    16060813
% Date:     Jun-08 2016, 1:15 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2016, Jochen Weber
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

% load default OLT
persistent ovl_defcolor;
if isempty(ovl_defcolor)
    ovl_defcolor = struct;
    try
        defolt = cell(1, 1);
        defolt{1} = xff([neuroelf_path('lut') '/Standard.olt']);
        defoltc = defolt{1}.C;
        ovl_defcolor.standard = defoltc.Colors;
        clearxffobjects(defolt);
    catch xfferror
        neuroelf_lasterr(xfferror);
        clearxffobjects(defolt);
        ovl_defcolor.standard = [255 * ones(1, 10), zeros(1, 10); ...
            [75:20:255, 75:20:255]; zeros(1, 10), 255:-20:75]';
    end
    ovl_defcolor.xstandard = [255 * ones(1, 40), zeros(1, 40); ...
        [60:5:255, 60:5:255]; zeros(1, 40), 248, 246, 243, 239, 235:-5:60]';
end

% argument check
if nargin < 1 || numel(xo) ~= 1 || ~xffisobject(xo, true, ...
    {'ava', 'cmp', 'fsmf', 'glm', 'hdr', 'head', 'mtc', 'smp', 'vmp', 'vtc'})
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
ft = lower(xo.S.Extensions{1});
hasmaps = (~any(strcmp(ft, {'ava', 'glm', 'hdr', 'head', 'mtc', 'vtc'})));

% depending on type
switch (ft)
    case 'ava'
        nmaps = numel(fieldnames(bc.Maps)) + (size(bc.Maps.CellMeans, ndims(bc.Maps.CellMeans)) - 1);
    case 'glm'
        nmaps = numel(bc.Predictor);
    case 'hdr'
        nmaps = size(bc.VoxelData, 4);
    case 'head'
        nmaps = numel(bc.Brick);
    case {'cmp', 'fsmf', 'smp', 'vmp'}
        nmaps = numel(bc.Map);
    case 'mtc'
        if ~isfield(bc.RunTimeVars, 'AvgMTC') || ~bc.RunTimeVars.AvgMTC
            nmaps = 1;
        else
            nmaps = numel(bc.RunTimeVars.ConditionNames);
        end
    case 'vtc'
        if ~isfield(bc.RunTimeVars, 'AvgVTC') || ~bc.RunTimeVars.AvgVTC
            nmaps = 1;
        else
            nmaps = numel(bc.RunTimeVars.ConditionNames);
        end
end
if ~hasmaps && (~isfield(bc.RunTimeVars, 'Map') || numel(bc.RunTimeVars.Map) ~= nmaps)
    bc.RunTimeVars.Map = repmat(struct( ...
        'Type', 15, 'LowerThreshold', 1, 'UpperThreshold', 5, 'Name', '', ...
        'RGBLowerThreshPos', [255, 0, 0], 'RGBUpperThreshPos', [255, 255, 0], ...
        'RGBLowerThreshNeg', [255, 0, 255], 'RGBUpperThreshNeg', [0, 0, 255], ...
        'UseRGBColor', 0, 'LUTName', '<default>', 'TransColorFactor', 1, ...
        'NrOfLags', 0, 'MinLag', 0, 'MaxLag', 0, 'CCOverlay', 0, ...
        'ClusterSize', 4, 'EnableClusterCheck', 0, 'UseValuesAboveThresh', 1, ...
        'DF1', 0, 'DF2', 0, 'ShowPositiveNegativeFlag', 3, ...
        'BonferroniValue', 0, 'NrOfFDRThresholds', 0, 'FDRThresholds', zeros(0, 3), ...
        'OverlayColors', []), 1, nmaps);
    if strcmp(ft, 'vtc') && ...
        isfield(bc.RunTimeVars, 'ConditionColors') && ...
        isequal(size(bc.RunTimeVars.ConditionColors), [nmaps, 16]) && ...
        isfield(bc.RunTimeVars, 'ConditionThresholds') && ...
       ~isempty(bc.RunTimeVars.ConditionThresholds) && ...
        size(bc.RunTimeVars.ConditionThresholds, 1) == nmaps && ...
        size(bc.RunTimeVars.ConditionThresholds, 3) == 2
        for cc = 1:nmaps
            bc.RunTimeVars.Map(cc).Type = 30;
            bc.RunTimeVars.Map(cc).LowerThreshold = bc.RunTimeVars.ConditionThresholds(cc, 1, 1);
            bc.RunTimeVars.Map(cc).UpperThreshold = bc.RunTimeVars.ConditionThresholds(cc, 1, 2);
            bc.RunTimeVars.Map(cc).Name = bc.RunTimeVars.ConditionNames{cc};
            bc.RunTimeVars.Map(cc).RGBLowerThreshPos = double(bc.RunTimeVars.ConditionColors(cc, 1:3));
            bc.RunTimeVars.Map(cc).RGBUpperThreshPos = double(bc.RunTimeVars.ConditionColors(cc, 5:7));
            bc.RunTimeVars.Map(cc).RGBLowerThreshNeg = double(bc.RunTimeVars.ConditionColors(cc, 9:11));
            bc.RunTimeVars.Map(cc).RGBUpperThreshNeg = double(bc.RunTimeVars.ConditionColors(cc, 13:15));
            bc.RunTimeVars.Map(cc).UseRGBColor = 1;
        end
    end
end
if nargin < 2 || isempty(maps)
    maps = 1:nmaps;
elseif ~isa(maps, 'double') || any(isinf(maps(:)) | isnan(maps(:)) | maps(:) < 1 | maps(:) > nmaps)
    error('neuroelf:xff:badArgument', 'Bad maps specification.');
else
    maps = round(real(maps(:)'));
end
if nargin < 3
    cspec = 'auto';
elseif numel(cspec) == 1 && xffisobject(cspec, true, 'olt')
    cspec = cspec.C.Colors;
elseif ischar(cspec)
    if ~any(strcmpi(cspec(:)', {'', 'auto', 'hsv', 'hsv4rgb', 'rgb', 'xauto', 'xhsv', 'xhsv4rgb'}))
        try
            cspecf = cell(1, 1);
            cspecf{1} = xff(cspec(:)');
            if ~xffisobject(cspecf{1}, true, 'olt')
                error('NO_OLT_FILE');
            end
            oltc = cspecf{1}.C;
            cspec = oltc.Colors;
            clearxffobjects(cspecf);
        catch xfferror
            clearxffobjects(cspecf);
            error('neuroelf:xff:badArgument', 'Invalid string cspec argument: %s.', xfferror.message);
        end
    else
        cspec = lower(cspec(:)');
    end
elseif isa(cspec, 'double')
    if ~isempty(cspec) && (size(cspec, 2) ~= 3 || mod(size(cspec, 1), 2) ~= 0 || ...
        size(cspec, 1) > 240 || any(isinf(cspec(:)) | isnan(cspec(:)) | cspec(:) < 0 | cspec(:) > 255))
        error('neuroelf:xff:badArgument', 'Invalid double cspec argument.');
    else
        cspec = round(real(cspec));
    end
    if isempty(cspec)
        cspec = [];
    end
else
    error('neuroelf:xff:badArgument', 'Invalid cspec argument.');
end

% iterate over selected maps
for mc = maps

    % get current color
    if hasmaps
        curcol = bc.Map(mc).OverlayColors;
    else
        curcol = bc.RunTimeVars.Map(mc).OverlayColors;
    end

    % double spec
    if ~ischar(cspec)

        % where to put
        if hasmaps
            bc.Map(mc).OverlayColors = cspec;
        else
            bc.RunTimeVars.Map(mc).OverlayColors = cspec;
        end
        continue;
    end

    % for auto and non-empty OverlayColors, do nothing
    if ~isempty(strfind(cspec, 'auto')) && ~isempty(curcol)
        continue;
    end

    % char spec with OLT support
    if any(strcmp(cspec, {'auto', 'hsv4rgb', 'xauto', 'xhsv4rgb'}))

        if (hasmaps && bc.Map(mc).UseRGBColor == 0) || ...
           (~hasmaps && bc.RunTimeVars.Map(mc).UseRGBColor == 0)
            if hasmaps
                lut = bc.Map(mc).LUTName(:)';
            else
                lut = bc.RunTimeVars.Map(mc).LUTName(:)';
            end
            if ~any(lut == '<' | lut == '>') && numel(lut) > 4
                try
                    lutobj = cell(1, 1);
                    lutobj{1} = xff(lut);
                    lutc = lutobj{1}.C;
                    if hasmaps
                        bc.Map(mc).OverlayColors = lutc.Colors;
                    else
                        bc.RunTimeVars.Map(mc).OverlayColors = lutc.Colors;
                    end
                    clearxffobjects(lutobj);
                    continue;
                catch xfferror
                    neuroelf_lasterr(xfferror);
                    clearxffobjects(lutobj);
                end
            end
            if hasmaps
                if cspec(1) == 'x'
                    bc.Map(mc).OverlayColors = ovl_defcolor.xstandard;
                else
                    bc.Map(mc).OverlayColors = ovl_defcolor.standard;
                end
            else
                if cspec(1) == 'x'
                    bc.RunTimeVars.Map(mc).OverlayColors = ovl_defcolor.xstandard;
                else
                    bc.RunTimeVars.Map(mc).OverlayColors = ovl_defcolor.standard;
                end
            end
            continue;
        end
    end

    % prepare 4x3 colors array
    try
        if hasmaps
            cols = [bc.Map(mc).RGBLowerThreshPos; bc.Map(mc).RGBUpperThreshPos; ...
                bc.Map(mc).RGBLowerThreshNeg; bc.Map(mc).RGBUpperThreshNeg];
        else
            cols = [bc.RunTimeVars.Map(mc).RGBLowerThreshPos; bc.RunTimeVars.Map(mc).RGBUpperThreshPos; ...
                bc.RunTimeVars.Map(mc).RGBLowerThreshNeg; bc.RunTimeVars.Map(mc).RGBUpperThreshNeg];
        end
        if ~isequal(size(cols), [4, 3]) || any(isinf(cols(:)) | isnan(cols(:)) | cols(:) < 0 | cols(:) > 255)
            error('BAD_COLORS');
        end
        cols = real(cols);
    catch xfferror
        neuroelf_lasterr(xfferror);
        cols = [255, 0, 0; 255, 255, 0; 255, 0, 255; 0, 0, 255];
    end

    % HSV?
    if ~isempty(strfind(cspec, 'hsv'))
        cols = rgb2hsv((1 / 255) .* cols);
        if abs(cols(1) - cols(2)) > 0.5
            cols(1) = cols(1) + sign(cols(2) - cols(1));
        end
        if abs(cols(3) - cols(4)) > 0.5
            cols(3) = cols(3) + sign(cols(4) - cols(3));
        end
    end

    % number of gradients
    if any(cspec == 'x')
        mm = repmat((0:1/39:1)', 1, 3);
    else
        mm = repmat((0:1/9:1)', 1, 3);
    end
    smm = size(mm, 1);

    % create array
    lut = [mm .* repmat(cols(2, :), smm, 1) + (1 - mm) .* repmat(cols(1, :), smm, 1); ...
        mm .* repmat(cols(4, :), smm, 1) + (1 - mm) .* repmat(cols(3, :), smm, 1)];

    % HSV ?
    if ~isempty(strfind(cspec, 'hsv'))
        lut(lut > 1) = lut(lut > 1) - 1;
        lut = 256 .* hsv2rgb(lut);
    end

    % store
    if hasmaps
        bc.Map(mc).OverlayColors = min(255, floor(lut + 0.001));
    else
        bc.RunTimeVars.Map(mc).OverlayColors = min(255, floor(lut + 0.001));
    end
end

% set back into object
xo.C = bc;
