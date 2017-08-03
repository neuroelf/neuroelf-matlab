function xo = vmp_ImportSPMMaps(xo, maps, opts)
% VMP::ImportSPMMaps  - import SPM maps from Analzye/NIftI files
%
% FORMAT:       vmp.ImportSPMMaps(maps [, opts])
%
% Input fields:
%
%       maps        filenames of spmT_* / spmF_* maps
%       opts        1x1 struct with optional fields
%        .interp    interpolation 'cubic', 'lanczos3', {'linear'}, 'nearest'
%        .maptype   1x1 or Mx1 map type, 'F', 'r', 't' (default from name)
%        .rpvmap    filename of RPV.hdr/nii filename (same space!)
%        .thresh    1x2 or Mx2 lower and upper thresholds (default: p=0.05)
%
% No output fields.
%
% Using: applyfdr, lsqueeze.

% Version:  v1.1
% Build:    16072114
% Date:     Jul-21 2016, 2:10 PM EST
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
lsqueeze = ne_methods.lsqueeze;

% check arguments
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'vmp') || isempty(maps) || ...
   (~ischar(maps) && ~xffisobject(maps, true, {'hdr', 'head'}) && ~iscell(maps))
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
res = bc.Resolution;
if isempty(bc.Map)
    error('neuroelf:xff:badArgument', 'The VMP must not be empty (for default settings).');
end
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ischar(maps)
    if any(size(maps) == 1)
        maps = {maps(:)'};
    else
        maps = cellstr(maps);
    end
elseif xffisobject(maps, true)
    maps = {maps};
end
if ~isfield(opts, 'interp') || ~ischar(opts.interp) || ...
   ~any(strcmpi(opts.interp(:)', {'cubic', 'lanczos3', 'linear', 'nearest'}))
    opts.interp = 'linear';
else
    opts.interp = lower(opts.interp(:)');
end
if ~isfield(opts, 'maptype') || ~ischar(opts.maptype) || ...
   (numel(opts.maptype) ~= 1 && numel(opts.maptype) ~= numel(maps))
    opts.maptype = 'a';
else
    opts.maptype = lower(opts.maptype(:)');
end
if numel(opts.maptype) == 1
    opts.maptype(2:numel(maps)) = opts.maptype(1);
end
if ~isfield(opts, 'rpvmap') || ~ischar(opts.rpvmap) || isempty(opts.rpvmap) || ...
    exist(opts.rpvmap(:)', 'file') ~= 2
    opts.rpvmap = '';
end
if ~isfield(opts, 'thresh') || (~isa(opts.thresh, 'double') && ~isa(opts.thresh, 'single')) || ...
    size(opts.thresh, 2) ~= 2 || (size(opts.thresh, 1) ~= 1 && size(opts.thresh, 1) ~= numel(maps)) || ...
    any(isinf(opts.thresh(:)) | isnan(opts.thresh(:)))
    opts.thresh = [0.05, 0.001];
else
    opts.thresh = abs(opts.thresh(:, 1:2, 1));
    opts.thresh = opts.thresh(:, 1:2);
end
if size(opts.thresh, 1) == 1
    opts.thresh = repmat(opts.thresh, [numel(maps), 1]);
end
for mc = numel(maps):-1:1
    if isempty(maps{mc}) || (~ischar(maps{mc}) && ...
        ~xffisobject(maps{mc}, true, {'hdr', 'head'})) || ...
       (ischar(maps{mc}) && exist(maps{mc}(:)', 'file') ~= 2)
        maps(mc) = [];
        opts.maptype(mc) = [];
        opts.thresh(mc, :) = [];
    end
    if ischar(maps{mc})
        maps{mc} = maps{mc}(:)';
        if numel(maps{mc}) > 4 && strcmpi(maps{mc}(end-3:end), '.img')
            if exist([maps{mc}(1:end-3) 'hdr'], 'file') == 2
                maps{mc} = [maps{mc}(1:end-3) 'hdr'];
            elseif exist([maps{mc}(1:end-3) 'HDR'], 'file') == 2
                maps{mc} = [maps{mc}(1:end-3) 'HDR'];
            else
                opts.maptype(mc) = [];
                opts.thresh(mc, :) = [];
                maps(mc) = [];
            end
        end
    end
end
if isempty(maps)
    error('neuroelf:xff:badArgument', 'Map file(s) not found.');
end
nobj = numel(maps);
nmaps = ones(nobj, 1);
hobj = false(1, nobj);
lobj = false(1, nobj);
mobj = cell(1, nobj);
rpvmap = false;
if isempty(opts.rpvmap)
    if xffisobject(maps{1}, true)
        opts.rpvmap = [fileparts(aft_FilenameOnDisk(maps{1})) '/RPV.hdr'];
        if exist(opts.rpvmap, 'file') ~= 2
            opts.rpvmap = [fileparts(aft_FilenameOnDisk(maps{1})) '/RPV.nii'];
        end
    else
        opts.rpvmap = [fileparts(maps{1}) '/RPV.hdr'];
        if exist(opts.rpvmap, 'file') ~= 2
            opts.rpvmap = [fileparts(maps{1}) '/RPV.nii'];
        end
    end
    if exist(opts.rpvmap, 'file') == 2
        try
            opts.rpvmap = xff(opts.rpvmap);
            if numel(opts.rpvmap) == 1 && xffisobject(opts.rpvmap, true, 'hdr')
                rpvmap = opts.rpvmap.C;
                if isempty(regexpi(rpvmap.DataHist.Description, '^spm_spm.*resels\s+per\s*voxel'))
                    error('BAD_RPV_CONTENT');
                end
                rpvmap = true;
            else
                error('BAD_RPV_HEAD_FILE');
            end
        catch xfferror
            neuroelf_lasterr(xfferror);
            if xffisobject(opts.rpvmap, true)
                clearxffobjects({opts.rpvmap});
            end
            opts.rpvmap = '';
            rpvmap = false;
        end
    end
end
for mc = 1:nobj
    try
        if ~xffisobject(maps{mc}, true)
            lobj(mc) = true;
            mobj{mc} = xff(maps{mc});
        else
            mobj{mc} = maps{mc};
        end
        if numel(mobj{mc}) ~= 1 || ~xffisobject(mobj{mc}, true, {'hdr', 'head'})
            error('BAD_HDR_HEAD_FILE');
        end
        tbc = mobj{mc}.C;
        if strcmpi(mobj{mc}.S.Extensions{1}, 'hdr')
            if isempty(regexpi(tbc.DataHist.Description, ...
                    '^spm (beta|contrast) - (\d+)\:\s*(.*)$')) && ...
                isempty(regexpi(tbc.DataHist.Description, ...
                    '^spm\{([a-z]+)_\[([0-9][0-9\.]+)(\,[0-9][0-9\.]+)?\]\}\s*\-?\s*(.*)$'))
                if opts.maptype(mc) == 'a'
                    [mapfp, mapff] = fileparts(maps{mc});
                    tbc.DataHist.Description = ['imported map: ' mapff];
                end
            end
            nmaps(mc) = size(tbc.VoxelData, 4);
            hobj(mc) = true;
        else
            nmaps(mc) = numel(tbc.Brick);
        end
    catch xfferror
        neuroelf_lasterr(xfferror);
        clearxffobjects(mobj(lobj));
        error('neuroelf:xff:badArgument', ...
            'Invalid map file (not an SPM result map): %s', maps{mc});
    end
end

% build coordinate indices for sampling
bb = struct('BBox', [bc.XStart, bc.YStart, bc.ZStart; bc.XEnd, bc.YEnd, bc.ZEnd], 'ResXYZ', res);
if res > 1 || bc.NativeResolutionFile
    bb.BBox(2, :) = bb.BBox(2, :) - 1;
else
    bb.BBox(2, :) = bb.BBox(2, :) + 0.5;
end

% add maps to structure once !
bc.Map = bc.Map(:);
onmaps = numel(bc.Map);
bc.Map(onmaps+1:onmaps+sum(nmaps)) = bc.Map(onmaps);

% FDR thresholds
fdrv = [0.1, 0.05, 0.04, 0.03, 0.02, 0.01, 0.005, 0.001]';

% RPV data
if rpvmap
    rpvimg = aft_SampleBVBox(opts.rpvmap, bb, 1, 'cubic');
    clearxffobjects({opts.rpvmap});
    rpvfwhm = ones(1, 3) .* (res ./ (mean(rpvimg(~isnan(rpvimg) & rpvimg > 1e-16 & rpvimg < 1)) .^ (1/3)));
    rpvdata = single(res ./ (rpvimg .^ (1/3)));
    rpvimg = single(rpvimg);
    rpvdata(isinf(rpvdata) | isnan(rpvdata) | rpvdata < 0) = 0;
    rpvdata(rpvdata > 0 & rpvdata < 1) = 1;
    rpvdata(rpvdata > (res * max(size(rpvdata)))) = res * max(size(rpvdata));
end

% iterate over map objects
tobj = onmaps;
for oc = 1:nobj

    % get content of map
    if isfield(mobj{oc}.H, 'GZIPfile') && ischar(mobj{oc}.H.GZIPfile) && ~isempty(mobj{oc}.H.GZIPfile)
        [sfn{1:2}] = fileparts(mobj{oc}.H.GZIPfile);
    else
        [sfn{1:2}] = fileparts(mobj{oc}.F);
    end
    if ~isempty(sfn{1})
        [pfn{1:2}] = fileparts(sfn{1});
        sfn{2} = [pfn{2} '/' sfn{2}];
    end
    cobj = mobj{oc}.C;
    cobo = cobj;
    if hobj(oc)
        if istransio(cobj.VoxelData)
            cobj.VoxelData = resolve(cobj.VoxelData);
        end
        cobj.VoxelData = single(cobj.VoxelData);
        cobj.VoxelData(isnan(cobj.VoxelData)) = single(0);
        bfv = sum(lsqueeze(any(cobj.VoxelData ~= 0, 4)));
        desc = cobj.DataHist.Description;
        [iti{1:3}] = regexpi(desc, '^spm (beta|contrast) - (\d+)\:\s*(.*)$');
        [dfi{1:3}] = regexpi(desc, ...
            '^spm\{([a-z]+)_\[([0-9][0-9\.]+)(\,[0-9][0-9\.]+)?\]\}\s*\-?\s*(.*)$');
        if ~any('bcfrtz' == opts.maptype(oc))
            if ~isempty(iti{3})
                opts.maptype(oc) = lower(desc(iti{3}{1}(1)));
                mapname = sprintf('%s - %s %s: %s', sfn{2}, ...
                    desc(iti{3}{1}(1, 1):iti{3}{1}(1, 2)), ...
                    desc(iti{3}{1}(2, 1):iti{3}{1}(2, 2)), ...
                    desc(iti{3}{1}(3, 1):iti{3}{1}(3, 2)));
            elseif ~isempty(dfi{3})
                opts.maptype(oc) = lower(desc(dfi{3}{1}(1)));
                mapname = sprintf('%s: %s', sfn{2}, ...
                    desc(dfi{3}{1}(4, 1):dfi{3}{1}(4, 2)));
            else
                mapname = sprintf('%s: %s', sfn{2}, desc);
            end
        else
            mapname = sprintf('%s: %s', sfn{2}, desc);
        end
        if ~any('bcfrtz' == opts.maptype(oc))
            opts.maptype(oc) = 'b';
        end
    else
        for omc = 1:nmaps(oc)
            if istransio(cobj.Brick(omc).Data)
                cobj.Brick(omc).Data = resolve(cobj.Brick(omc).Data);
            end
            cobj.Brick(omc).Data = single(cobj.Brick(omc).Data);
            cobj.Brick(omc).Data(isnan(cobj.Brick(omc).Data)) = single(0);
        end
        bfv = sum(lsqueeze(any(cat(4, cobj.Brick.Data) ~= 0, 4)));
    end
    mobj{oc}.C = cobj;

    % parse individual maps
    for omc = 1:nmaps(oc)

        % get target position
        tobj = tobj + 1;

        % set DF parameters
        if ~hobj(oc)
            if ~isempty(cobj.Brick(omc).FuncParams)
                bc.Map(tobj).DF1 = cobj.Brick(omc).FuncParams(1);
                if numel(cobj.Brick(omc).FuncParams) > 1
                    bc.Map(tobj).DF2 = cobj.Brick(omc).FuncParams(2);
                end
            end
            mapname = [sfn{2} ' - ' cobj.Brick(omc).Label];
        end

        % fill in map values
        bc.Map(tobj).DF2 = 0;
        switch (opts.maptype(oc))
            case {'b', 'c'}
                bc.Map(tobj).Type = 11;
            case 'f'
                bc.Map(tobj).Type = 4;
                if hobj(oc) && ~isempty(dfi) && ~isempty(dfi{3}) && ...
                    dfi{3}{1}(3, 2) >= dfi{3}{1}(3, 1)
                    bc.Map(tobj).DF2 = str2double(desc(dfi{3}{1}(3, 1):dfi{3}{1}(3, 2)));
                else
                    bc.Map(tobj).DF2 = 1;
                end
            case 'r'
                bc.Map(tobj).Type = 2;
            case 't'
                bc.Map(tobj).Type = 1;
            case 'z'
                bc.Map(tobj).Type = 12;

            % for HEAD files
            otherwise
                if cobj.Brick(omc).FuncType == 0
                    bc.Map(tobj).Type = 11;
                else
                    bc.Map(tobj).Type = cobj.Brick(omc).FuncType;
                end
        end
        if hobj(oc)
            if any('frt' == opts.maptype(oc))
                if ~isempty(dfi{3})
                    bc.Map(tobj).DF1 = str2double(desc(dfi{3}{1}(2, 1):dfi{3}{1}(2, 2)));
                else
                    bc.Map(tobj).DF1 = 120;
                end
                bc.Map(tobj).FDRThresholds = [fdrv, ...
                    ne_methods.applyfdr(double(cobj.VoxelData(:)), opts.maptype(oc), fdrv, ...
                    bc.Map(tobj).DF1, bc.Map(tobj).DF2, true)];
            elseif (opts.maptype(oc) == 'z');
                bc.Map(tobj).DF1 = 5e3;
                bc.Map(tobj).FDRThresholds = [fdrv, ...
                    ne_methods.applyfdr(double(cobj.VoxelData(:)), 't', fdrv, ...
                    bc.Map(tobj).DF1, bc.Map(tobj).DF2, true)];
            else
                bc.Map(tobj).DF1 = 1;
                bc.Map(tobj).FDRThresholds = [1, 1e9, 1e9];
            end
        else
            bc.Map(tobj).FDRThresholds = zeros(0, 3);
        end
        bc.Map(tobj).NrOfFDRThresholds = size(bc.Map(tobj).FDRThresholds, 1);
        bc.Map(tobj).Name = mapname;
        bc.Map(tobj).BonferroniValue = bfv;
        bc.Map(tobj).UnknownValue = -1;
        bc.Map(tobj).VMPData = single(aft_SampleBVBox(mobj{oc}, bb, omc, opts.interp));
        bc.Map(tobj).VMPDataCT = [];

        % also store information from RPV?
        if rpvmap
            bc.Map(tobj).RunTimeVars.FWHMResEst = rpvfwhm;
            bc.Map(tobj).RunTimeVars.FWHMResImg = rpvdata;
            bc.Map(tobj).RunTimeVars.RPVImg = rpvimg;
        end
    end

    % reduce memory footprint (revert to transio)
    mobj{oc}.C = cobo;
end

% clear temporary objects
clearxffobjects(mobj(lobj));

% put content in new object
bc.NrOfMaps = numel(bc.Map);
xo.C = bc;
