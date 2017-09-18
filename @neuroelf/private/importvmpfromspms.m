function vmp = importvmpfromspms(maps, maptype, bbox, res, imeth)
% importvtcfromanalyze  - import several Analzye SPM map files
%
% FORMAT:       vmp = importvmpfromspms([maps [, type [, bbox, res [, imeth]]]])
%
% Input fields:
%
%       maps        char or cell array with SPM map filenames
%       type        type of map ('a' for auto, default)
%       bbox        optional 2x3 bounding box (default: small TAL box)
%       res         optional resolution (default: 3)
%       imeth       interpolation 'cubic', 'lanczos3', {'linear'}, 'nearest'
%
% Output fields:
%
%       vmp         created VMP object

% Version:  v1.1
% Build:    17091810
% Date:     Sep-18 2017, 10:34 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2017, Jochen Weber
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

% argument check
if nargin < 1 || ...
    isempty(maps) || ...
   (~iscell(maps) && ...
    ~ischar(maps))
    if mainver > 6
        msargs = {'MultiSelect', 'on'};
    else
        msargs = {};
    end
    vmp = [];
    [mapf, mapp] = uigetfile( ...
        {'*.hdr;*.nii;*.nii.gz', 'SPM result map files (*.hdr, *.nii, *.nii.gz)'; ...
         '*.mat', 'SPM result file (*.mat)'; ...
         '*.head', 'AFNI head files (*.head)'}, ...
        'Please select the SPM map files or SPM.mat file to import...', '*.hdr', msargs{:});
    if isequal(mapf, 0) || ...
        isequal(mapp, 0)
        return;
    end
    if ~iscell(mapf)
        maps = {strrep([strrep(mapp, '\', '/') '/' mapf], '//', '/')};
    else
        maps = mapf;
        for mc = 1:numel(maps)
            maps{mc} = strrep([strrep(mapp, '\', '/') '/' maps{mc}], '//', '/');
        end
    end
    mc = 1;
    maps = maps(:);
    while mc <= numel(maps)
        if ~isempty(regexpi(maps{mc}, '\.mat$'))
            try
                spmc = load(maps{mc});
                if ~isfield(spmc, 'SPM')
                    error('NOT_AN_SPM_FILE');
                end
                spmp = fileparts(maps{mc});
                coni = {spmc.SPM.xCon.Vspm};
                for smc = 1:numel(coni)
                    [conip, conif, conie] = fileparts(coni{smc}.fname);
                    coni{smc} = [spmp '/' conif conie];
                end
                maps = [maps(1:mc-1); coni(:); maps(mc+2:end)];
                mc = mc + numel(coni) - 1;
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
                return;
            end
        end
        mc = mc + 1;
    end
end
if ~iscell(maps)
    if size(maps, 1) > 1
        maps = cellstr(maps);
    else
        maps = {maps};
    end
end
if nargin < 2 || ...
   ~ischar(maptype)
    maptype = 'a';
end
if nargin < 3 || ...
   ~isa(bbox, 'double') || ...
   ~isequal([2, 3], size(bbox))
    bbox = [44, 38, 44; 241, 193, 211];
end
if nargin < 4 || ...
   ~isa(res, 'double') || ...
    numel(res) ~= 1 || ...
   ~any((1:12) == res)
    res = 3;
end
if nargin < 5 || ...
   ~ischar(imeth) || ...
   ~any(strcmpi(imeth(:)', {'cubic', 'lanczos3', 'linear', 'nearest'}))
    imeth = 'linear';
else
    imeth = lower(imeth(:)');
end

% error handling
try
    vmps = cell(1, 1);
    vmp = newnatresvmp(bbox, res);
    vmps{1} = vmp;
    opts = struct;
    opts.interp = imeth;
    opts.maptype = maptype;
    vmp.ImportSPMMaps(maps, opts);
    vmp.Map = vmp.Map(2:end);
    vmp.NrOfMaps = numel(vmp.Map);
catch ne_eo;
    clearxffobjects(vmps);
    rethrow(ne_eo);
end
