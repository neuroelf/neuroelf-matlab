function nlf = aft_ConvertToNLF(xo, opts)
% AFT::ConvertToNLF  - convert an object to NLF (if appropriate)
%
% FORMAT:       obj.ConvertToNLF([opts]);
%
% Input fields:
%
%       opts        optional settings
%        .nodata    flag, if data is transio, use reference (default: false)
%        .vmr16     flag, choose VMRData16 (if available, default: false)
%
% No output fields.
%
% TYPES: GLM, SMP, SRF, VMP, VMR, VTC
%
% Using: mltype.

% Version:  v1.1
% Build:    16020213
% Date:     Feb-02 2016, 1:59 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2014, 2016, Jochen Weber
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

% check arguments
if numel(xo) ~= 1 || ~xffisobject(xo, true, {'glm', 'smp', 'srf', 'vmp', 'vmr', 'vtc'})
    error('neuroelf:xff:badArguments', ...
        'This type of object does not support To-NLF-conversion.');
end

% options
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'nodata') || ~islogical(opts.nodata) || numel(opts.nodata) ~= 1
    opts.nodata = false;
end
if ~isfield(opts, 'vmr16') || ~islogical(opts.vmr16) || numel(opts.vmr16) ~= 1
    opts.vmr16 = false;
end

% get content and object type
bc = xo.C;
ext = lower(xo.S.Extensions{1});

% get filename particles
[hfdir, hffile, hfext] = fileparts(xo.F);
if isempty(hffile)
    hffile = 'unsaved';
end
if isempty(hfext)
    hfext = sprintf('.%s', ext);
end

% create and get NLF content, then destroy object
nlf = xff('new:nlf');
nc = nlf.C;

% depending on extension
switch (ext)

    % for VMRs
    case 'vmr'

        % get bounding box
        bbox = aft_BoundingBox(xo);

        % ensure intent code
        nc.Intent = 'a3d';
        nc.NrOfDims = 3;
        nc.DimMeaning = 'xyz.....';
        nc.DimUnit = 'mmm.....';

        % put data into appropriate fields
        if opts.vmr16 && ...
           ~isempty(bc.VMRData16)
            nc.Data = bc.VMRData16;
            nc.DataType = ne_methods.mltype('uint16');
        else
            nc.Data = bc.VMRData;
            nc.DataType = ne_methods.mltype('uint8');
            nc.A3D.NrOfSVCs = 30;
            nc.A3D.SVColors = [226:255; ...
                 255 .* ones(1, 10), zeros(1, 10), 255 .* ones(1, 10); ...
                 75:20:255, 75:20:255, 255 .* ones(1, 10); ...
                 zeros(1, 10), 255:-20:75, 255 .* ones(1, 10)]';
        end

        % some more settings
        nc.Size = [size(nc.Data), ones(1, 5)];
        nc.Resolution = [bc.VoxResX, bc.VoxResY, bc.VoxResZ, ones(1, 5)];
        nc.ScalingSlope = 1;
        nc.ScalingIntercept = 0;
        nc.GlobalTransform = bbox.QuatB2T;
        nc.Name = [hffile, hfext];
        nc.A3D.NrOfTransforms = numel(bc.Trf);
        for tc = 1:numel(bc.Trf)
            nc.A3D.Transform(tc).Type = 32 + bc.Trf(tc).TypeOfSpatialTransformation;
            nc.A3D.Transform(tc).NrOfValues = numel(bc.Trf(tc).TransformationValues);
            nc.A3D.Transform(tc).Values = bc.Trf(tc).TransformationValues(:);
        end
end

% for transio data
if istransio(nc.Data)

    % keep reference data
    if opts.nodata

        % update filename, etc.
        nc.DataType = ne_methods.mltype(datatype(nc.Data));
        nc.DataExternal = 1;
        nc.DataExternalEndian = littleendian(nc.Data);
        nc.DataFile = filename(nc.Data);
        nc.DataOffset = offset(nc.Data);

    % resolve data
    else
        nc.Data = resolve(nc.Data);
    end
end

% put into global position of current xo, replacing original object
nc.SourceFilename = xo.F;
nc.RunTimeVars = bc.RunTimeVars;

% set content and keep handles (borderline)
nlf.C = nc;
nlf.H = xo.H;
