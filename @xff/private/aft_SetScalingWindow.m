function xo = aft_SetScalingWindow(xo, win, uhist)
% AFT::SetScalingWindow  - set the ScalingWindow (in RunTimeVars)
%
% FORMAT:       [obj = ] obj.SetScalingWindow([win [, uhist]])
%
% Input fields:
%
%       win         optional 1x2 min/max values for visualization,
%                   if not given or invalid, define values from data/type
%       uhist       force update scaling histogram (false)
%
% TYPES: AMR, DMR, FMR, HDR, HEAD, MGH, MSK, VMR, VTC
%
% Using: computescalingwindow.

% Version:  v1.1
% Build:    16031614
% Date:     Mar-16 2016, 2:35 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

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

% neuroelf library
global ne_methods;

% argument check
if nargin < 1 || numel(xo) ~= 1 || ~xffisobject(xo, true, {'amr', 'dmr', 'fmr', 'hdr', 'head', 'mgh', 'msk', 'vmr', 'vtc'})
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
if nargin < 2 || ~isa(win, 'double') || numel(win) ~= 2 || any(isnan(win))
    win = [-Inf, Inf];
end
if nargin < 3 || ~islogical(uhist) || numel(uhist) ~= 1
    uhist = false;
end

% nothing to do
if ~any(isinf(win)) && ~uhist

    % get content
    bc = xo.C;

    % limits not yet established
    if ~isfield(bc.RunTimeVars, 'ScalingWindowLim')
        aft_SetScalingWindow(xo);
        bc = xo.C;
    end

    % update window with limits, set and return
    bc.RunTimeVars.ScalingWindow = ...
       [max(bc.RunTimeVars.ScalingWindowLim(1), win(1)), ...
        min(bc.RunTimeVars.ScalingWindowLim(2), win(2))];

    xo.C = bc;
    return;
end

% otherwise, the extension is needed (to access the correct data field)
bc = xo.C;
ft = lower(xo.S.Extensions{1});
sci = 0;
scs = 1;

% depending on type
d16 = [];
switch (ft)

    % AMR -> AMRData
    case 'amr'

        % data in .Slice(s).AMRData
        d = bc.Slice(1).AMRData(:, :);
        d(1, 1, numel(bc.Slice)) = d(1);
        for sc = 2:numel(bc.Slice)
            d(:, :, sc) = bc.Slice(sc).AMRData(:, :);
        end

    % DMR -> DWIData
    case 'dmr'

        % data in .DWIData
        d = bc.DWIData(:);

    % FMR -> STCData
    case 'fmr'

        % data in .Slice(s).STCData
        d = bc.Slice(1).STCData(:, :, :, :);
        d(1, 1, 1, numel(bc.Slice)) = d(1, 1, 1, end);
        for sc = 2:numel(bc.Slice)
            d(:, :, :, sc) = bc.Slice(sc).STCData(:, :, :);
        end

    % HDR -> VoxelData
    case 'hdr'

        % special datatypes
        if any(bc.ImgDim.DataType == [128, 2304])
            d = bc.VoxelDataRGBA(:, :, :, :, 1:3);
        else
            % data in .VoxelData
            d = bc.VoxelData(:);
        end
        if any([2, 4, 8, 130, 132, 136, 256, 512, 768] == bc.ImgDim.DataType)
            if bc.ImgDim.ScalingIntercept ~= 0
                sci = bc.ImgDim.ScalingIntercept;
            end
            if all([0, 1] ~= bc.ImgDim.ScalingSlope)
                scs = bc.ImgDim.ScalingSlope;
            end
        end

    % HEAD -> Data
    case 'head'

        % data in .Brick(b).Data
        d = bc.Brick(1).Data(:, :, :);
        ds = size(d);
        d(1, 1, 1, numel(bc.Brick)) = d(1);
        scfs = ones(numel(bc.Brick), 1);
        if all([0, 1] ~= bc.Brick(1).ScalingFactor)
            scfs(1) = bc.Brick(1).ScalingFactor;
        end
        for sc = 2:numel(bc.Brick)
            d(:, :, :, sc) = reshape(bc.Brick(sc).Data(:), ds);
            if all([0, 1] ~= bc.Brick(sc).ScalingFactor)
                scfs(sc) = bc.Brick(sc).ScalingFactor;
            end
        end
        if all(scfs == scfs(1))
            scs = scfs(1);
        else
            d = reshape(reshape(d, prod(ds), numel(scfs)) * ...
                sparse(diag(scfs)), [ds, numel(scfs)]);
        end

    % MGH -> MGHData
    case 'mgh'

        % data in .VMRData
        d = bc.MGHData(:);

    % MSK -> special treatment (binary data!)
    case 'msk'
        msk = sum(bc.Mask(:) == 1);
        win = [0, 4 / 3];
        mmh = [zeros(8, 1); numel(bc.Mask) - msk; zeros(183, 1); msk; zeros(63, 1)];
        bc.RunTimeVars.ScalingWindow = win;
        bc.RunTimeVars.ScalingWindowLim = bc.RunTimeVars.ScalingWindow;
        bc.RunTimeVars.ScalingHist = mmh;
        xo.C = bc;
        return;

    % VMR -> VMRData
    case 'vmr'

        % data in .VMRData
        d = bc.VMRData(:);
        if isequal(size(bc.VMRData), size(bc.VMRData16))
            d16 = bc.VMRData16(:);
        end

    % VTC -> VTCData
    case 'vtc'

        % data in .VTCData
        d = bc.VTCData(:);
end

% replace Inf/NaN values
if isa(d, 'single') || isa(d, 'double')
    if size(d, 4) > 1
        ds3 = size(d, 1) * size(d, 2) * size(d, 3);
        for fdc = 1:size(d, 4)
            dp = d(:, :, :, fdc);
            dp = find(isinf(dp(:)) | isnan(dp(:)));
            if ~isempty(dp)
                d(ds3 * (fdc - 1) + dp) = meannoinfnan(dp(:));
            end
        end
    else
        mmd = ~isinf(d);
        mmd = mmd & ~isnan(d);
        if ~all(mmd(:))
            mmm = mean(d(mmd));
            if isnan(mmm)
                mmm = 0;
            end
            mmd = ~mmd;
            d(mmd) = mmm;
        end
    end
end

% sub-function
[win, mmh] = ne_methods.computescalingwindow(d, win, sci, scs);

% set back into object
bc.RunTimeVars.ScalingWindow = win;
bc.RunTimeVars.ScalingWindowLim = bc.RunTimeVars.ScalingWindow;
bc.RunTimeVars.ScalingHist = mmh;

% also compute for V16
if ~isempty(d16)
    [win16, mmh16] = ne_methods.computescalingwindow(d16, [0, Inf], 0, 1);
    bc.RunTimeVars.ScalingWindow16 = win16;
    bc.RunTimeVars.ScalingWindowLim16 = bc.RunTimeVars.ScalingWindow16;
    bc.RunTimeVars.ScalingHist16 = mmh16;
end

% for RGBA HDR, set Lim to [0, 255]
if strcmp(ft, 'hdr') && any(bc.ImgDim.DataType == [128, 2304])
    bc.RunTimeVars.ScalingWindowLim = [0, 255];
end

% set back to global storage
xo.C = bc;
