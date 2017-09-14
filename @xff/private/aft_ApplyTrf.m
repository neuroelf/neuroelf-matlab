function xo2 = aft_ApplyTrf(xo, trf, opts)
% AFT::ApplyTrf  - reslice data with a set of transformations
%
% FORMAT:       trfobj = obj.ApplyTrf(trfs [, opts]);
%
% Input fields:
%
%       trfs        single TRF object or cell of TRF (+TAL) objects
%       opts        optional settings
%        .bbox      2x3 alternative bounding box (default: same)
%        .method    resampling, one of {'cubic'}, 'lanczos3', 'linear', 'nearest'
%        .res       1x1 alternative resolution (default: same)
%
% Output fields:
%
%       trfobj      transformed object
%
% TYPES: CMP, GLM, MSK, VMP, VTC
%
% Using: applyspmsnc, bvcoordconv, findfirst, flexinterpn_method,
%        indexarraynb, limitrangec, lsqueeze, samplefmrspace.

% neuroelf library
global ne_methods;
acpc2tal           = ne_methods.acpc2tal;
flexinterpn_method = ne_methods.flexinterpn_method;

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, {'cmp', 'glm', 'msk', 'vmp', 'vtc'}) || ...
    isempty(trf) || (~xffisobject(trf, true) && ~iscell(trf))
    error('neuroelf:xff:badArgument', 'Bad or missing argument.');
end
if ~iscell(trf)
    trf = {trf};
else
    trf = trf(:);
end
for tc = 1:numel(trf)
    if ~xffisobject(trf{tc}, true, {'tal', 'trf'})
        error('neuroelf:xff:badArgument', 'Bad or missing argument.');
    end
end

% check object for type
xot = aft_Filetype(xo);
if strcmpi(xot, 'glm') && xo.C.ProjectType ~= 1
    error('neuroelf:xff:badArgument', 'Invalid GLM object.');
end

% get bounding box from object
bbox = aft_BoundingBox(xo);
itrf = bbox.QuatT2B;

% options
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'bbox') || ~isa(opts.bbox, 'double') || ~isequal(size(opts.bbox), [2, 3]) || ...
    any(isinf(opts.bbox(:)) | isnan(opts.bbox(:)))
    opts.bbox = bbox.BBox;
end
bb = opts.bbox;
if ~isfield(opts, 'method') || ~ischar(opts.method) || isempty(opts.method) || ...
   ~any(strcmpi(opts.method(:)', {'cubic', 'lanczos3', 'linear', 'nearest'}))
    opts.method = 'cubic';
else
    opts.method = lower(opts.method(:)');
end
if ~isfield(opts, 'res') || ~isa(opts.res, 'double') || numel(opts.res) ~= 1 || ...
    isinf(opts.res) || isnan(opts.res) || opts.res < 1
    opts.res = bbox.ResXYZ(1);
else
    opts.res = floor(opts.res);
end
r = opts.res;

% copy object
xo2 = aft_CopyObject(xo);

% generate sampling grid (values at which to sample)
[xg, yg, zg] = ndgrid(bb(1, 3):r:bb(2, 3)-.5, bb(1, 1):r:bb(2, 1)-.5, bb(1, 2):r:bb(2, 2)-.5);
nx = size(xg, 1);
ny = size(xg, 2);
nz = size(xg, 3);

% pack together in required system
xyz = [xg(:), yg(:), zg(:)];

% adapt bounding box and resolution
xo2.C.XStart = yg(1);
xo2.C.YStart = zg(1);
xo2.C.ZStart = xg(1);
xo2.C.XEnd = yg(end) + r;
xo2.C.YEnd = zg(end) + r;
xo2.C.ZEnd = xg(end) + r;

% apply transformations in backwards order
for tc = numel(trf):-1:1
    t = trf{tc};

    % TAL transform
    if isfield(t.C, 'AC')
        xyz = acpc2tal(xyz, t, true);
        continue;
    end
    
    % regular transformation
    xyz = [128 - xyz(:, [2, 3, 1]), ones(size(xyz, 1), 1)];
    xyz = xyz * inv(t.C.TFMatrix)';
    xyz = 128 - xyz(:, [3, 1, 2]);
end

% depending on object
switch lower(xot)
    
    % CMP/VMP
    case {'cmp', 'vmp'}
        if lower(xot(1)) == 'c'
            mapfield = 'CMPData';
        else
            mapfield = 'VMPData';
            xo2.C.NativeResolutionFile = 1;
            xo2.C.FileVersion = 6;
        end

        % iterate over maps
        for mc = 1:numel(xo2.C.Map)
            resmap = aft_SampleData3D(xo, 128 - xyz, struct('mapvol', mc, 'method', opts.method));
            xo2.C.Map(mc).(mapfield) = single(permute(reshape(resmap, [nx, ny, nz]), [2, 3, 1]));
            xo2.C.Map(mc).([mapfield 'CT']) = [];
        end

    % GLM
    case {'glm'}

        % RFX
        if xo2.C.ProjectTypeRFX > 0
            nm = size(xo2.C.GLMData.Subject(1).BetaMaps, 4);
            for sc = 1:numel(xo2.C.GLMData.Subject)
                xo2.C.GLMData.Subject.BetaMaps = single(zeros([ny, nz, nx, nm]));
            end
            tnm = aft_NrOfMaps(xo);
            sc = 1;
            smc = 1;
            for mc = 1:tnm
                resmap = aft_SampleData3D(xo, 128 - xyz, struct('mapvol', mc, 'method', opts.method));
                xo2.C.GLMData.Subject(sc).BetaMaps(:, :, :, smc) = ...
                    single(permute(reshape(resmap, [nx, ny, nz]), [2, 3, 1]));
                smc = smc + 1;
                if smc > nm
                    smc = 1;
                    sc = sc + 1;
                end
            end
            resmap = aft_SampleData3D(xo, 128 - xyz, struct('mapvol', tnm + 1, 'method', opts.method));
            xo2.C.GLMData.RFXGlobalMap = single(permute(reshape(resmap, [nx, ny, nz]), [2, 3, 1]));

        % FFX
        else
            emap = single(zeros([ny, nz, nx]));
            tnm = size(xo.C.GLMData.BetaMaps, 4);
            xo2.C.GLMData.MultipleRegressionR = emap;
            xo2.C.GLMData.MCorrSS = emap;
            xo2.C.GLMData.BetaMaps = repmat(emap, [1, 1, 1, tnm]);
            xo2.C.GLMData.XY = repmat(emap, [1, 1, 1, tnm]);
        end

    % MSK
    case {'msk'}
        if opts.method(1) ~= 'l' && opts.method(1) ~= 'n'
            opts.method = 'linear';
        end
        resmap = aft_SampleData3D(xo, 128 - xyz, struct('mapvol', 1, 'method', opts.method));
        xo2.C.MSKData = uint8(permute(reshape(resmap, [nx, ny, nz]), [2, 3, 1]) >= 0.5);

    % VTC
    case {'vtc'}
        nvol = size(xo.C.VTCData, 1);
        xo2.C.VTCData = xo.C.VTCData(1);
        xo2.C.VTCData(:) = 0;
        xo2.C.VTCData = repmat(xo2.C.VTCData, [nvol, ny, nz, nx]);
        for vc = 1:nvol
            resmap = aft_SampleData3D(xo, 128 - xyz, struct('mapvol', vc, 'method', opts.method));
            xo2.C.VTCData(vc, :, :, :) = permute(reshape(resmap, [1, nx, ny, nz]), [1, 3, 4, 2]);
        end
end
