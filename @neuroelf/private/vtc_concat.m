function [varargout] = vtc_concat(targetfile, vtcs, varargin)
% vtc_concat  - concatenate VTCs
%
% FORMAT:       newvtc = vtc_concat(targetfile, vtclist [, options]);
%         or    newvtc = vtc_concat(targetfile, vtc1, vtc2, ... [, options]);
%
% Input fields:
%
%       targetfile  filename of VTC to write
%       vtclist     cell array with VTCs to concatenated
%       vtc1, ...   single VTCs to concatenated
%       options     1x1 struct with optional settings
%        .datatype  target VTC DataType (default: keep, if filter/trans: 2)
%        .filtsec   filter cut-off in seconds (default: Inf)
%        .filttype  filter type, either of {'dct'}, 'fourier', or 'poly'
%        .prt       filename of new protocol file, otherwise keep first
%        .prts      cell array with PRTs to concatenate
%        .trans     transformation to apply on each VTC on reading
%                   either of 'psc' or 'z', leave unset for none
%        .volsel    cell array with volumes to keep (after filtering)
%        .xconfound cell array of files with regressors to remove
%
% Output fields:
%
%       newvtc      object handle to newly written VTC

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:30 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2012, 2014, 2016, Jochen Weber
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
if nargin < 2 || ...
   ~ischar(targetfile) || ...
    numel(targetfile) < 5 || ...
   (~all(isxff(vtcs, 'vtc')) && ...
    ~iscell(vtcs)) || ...
    isempty(vtcs)
    error( ...
        'neuroelf:BadArgument', ...
        'You must give a valid target filename and some VTCs.' ...
    );
end
if all(isxff(vtcs, 'vtc')) && ...
    numel(vtcs) > 1
    vtcc = cell(1, numel(vtcs));
    for vc = 1:numel(vtcs)
        vtcc{vc} = vtcs(vc);
    end
    vtcs = vtcc;
end
if ~iscell(vtcs)
    vtcs = {vtcs};
end
for n = 3:nargin
    if isxff(varargin{n - 2}, 'vtc')
        vtcs{end + 1} = varargin{n - 2};
    end
end
nv = numel(vtcs);
if nv == 1
    error( ...
        'neuroelf:TooFewArguments', ...
        'You need at least two VTCs to concatenate.' ...
    );
end

% options ?
datatype = 0;
filtsec = Inf;
filttype = 'dct';
newprt = '';
prts = {};
transopt = 0;
volsel = {};
xconfound = {};
if nargin > 2 && ...
    isstruct(varargin{end}) && ...
    numel(varargin{end}) == 1
    options = varargin{end};
else
    options = struct;
end
if isfield(options, 'datatype') && ...
    isa(options.datatype, 'double') && ...
   ~isinf(options.datatype) && ...
   ~isnan(options.datatype) && ...
    any([1, 2] == options.datatype)
    datatype = options.datatype;
end
if isfield(options, 'filtsec') && ...
    isa(options.filtsec, 'double') && ...
    numel(options.filtsec) == 1 && ...
   ~isnan(options.filtsec) && ...
    options.filtsec > 30
    filtsec = options.filtsec;
end
if isfield(options, 'filttype') && ...
    ischar(options.filttype) && ...
    any(strcmpi(options.filttype(:)', {'f', 'fourier'}))
    filttype = 'fourier';
elseif isfield(options, 'filttype') && ...
    ischar(options.filttype) && ...
    any(strcmpi(options.filttype(:)', {'p', 'poly'}))
    filttype = 'poly';
end
if isfield(options, 'prt') && ...
    ischar(options.prt)
    newprt = options.prt(:)';
end
if isfield(options, 'prts') && ...
    iscell(options.prts) && ...
    numel(options.prts) == nv && ...
    all(cellfun(@ischar, options.prts(:))) && ...
   ~any(cellfun('isempty', options.prts(:)))
    prts = options.prts(:);
end
for vc = 1:numel(prts)
    try
        prts{vc} = xff(prts{vc});
        if ~isxff(prts{vc}, 'prt')
            error('neuroelf:BadArgument', ...
                'Not a valid PRT for run %d.', ...
                vc ...
            );
        end
    catch ne_eo;
        clearxffobjects(prts);
        rethrow(ne_eo);
    end
end
if isfield(options, 'trans') && ...
    ischar(options.trans) && ...
   ~isempty(options.trans)
    switch (lower(options.trans(:)'))
        case {'psc', '%'}
            transopt = 1;
        case {'z'}
            transopt = 2;
    end
end
if isfield(options, 'volsel') && ...
    isa(options.volsel, 'double') && ...
   ~any(isinf(options.volsel(:)) | isnan(options.volsel(:)) | options.volsel(:) < 1)
    volsel = repmat({unique(round(options.volsel(:)))}, nv, 1);
elseif isfield(options, 'volsel') && ...
    iscell(options.volsel) && ...
    numel(options.volsel) == nv && ...
    all(cellfun(@isnumeric, options.volsel(:))) && ...
   ~any(cellfun('isempty', options.volsel(:)))
    volsel = options.volsel(:);
end
if isfield(options, 'xconfound') && ...
    iscell(options.xconfound) && ...
    numel(options.xconfound) == nv && ...
   ~any(cellfun('isempty', options.xconfound(:)))
    xconfound = options.xconfound(:);
end
for vc = 1:numel(xconfound)
    if ischar(xconfound{vc})
        try
            if ~isempty(regexpi(xconfound{vc}(:)', '\.mat$'))
                xcm = load(xconfound{vc}(:)');
                xcmf = fieldnames(xcm);
                if numel(xcmf) ~= 1
                    error( ...
                        'neuroelf:BadArgument', ...
                        'Invalid MAT-file for xconfound of run %d.', ...
                        vc ...
                    );
                end
                xcm = xcm.(xcmf{1});
            elseif ~isempty(regexpi(xconfound{vc}(:)', '\.(rtc|sdm)$'))
                xcmf = xff(xconfound{vc}(:)');
                if ~isxff(xcmf, 'rtc')
                    if isxff(xcmf, true)
                        xcmf.ClearObject;
                    end
                    error( ...
                        'neuroelf:BadArgument', ...
                        'Not an RTC/SDM file for run %d.', ...
                        vc ...
                    );
                end
                xcm = xcmf.SDMMatrix;
                xcmf.ClearObject;
            elseif ~isempty(regexpi(xconfound{vc}(:)', '\.txt$'))
                xcm = load(xconfound{vc}(:)');
            else
                error( ...
                    'neuroelf:BadArgument', ...
                    'Unsupported xconfound filename/type for run %d.', ...
                    vc ...
                );
            end
        catch ne_eo;
            clearxffobjects(prts);
            rethrow(ne_eo);
        end
        if isempty(xcm) || ...
            any(isinf(xcm(:)) | isnan(xcm(:)))
            clearxffobjects(prts);
            error( ...
                'neuroelf:BadArgument', ...
                'Invalid numeric content of xconfound for run %d.', ...
                vc ...
            );
        end
        xconfound{vc} = xcm;
    elseif ~isa(xconfound{vc}, 'double') || ...
        any(isinf(xconfound{vc}(:)) | isnan(xconfound{vc}(:)))
        clearxffobjects(prts);
        error( ...
            'neuroelf:BadArgument', ...
            'Invalid xconfound for run %d.', ...
            vc ...
        );
    end
end
if datatype == 0 && ...
   (~isinf(filtsec) || ...
    transopt > 0)
    datatype = 2;
elseif datatype == 0
    datatype = vtcs{1}.DataType;
end

% check vtcs argument for filenames
loadvtcs = false(1, nv);
if isempty(prts)
    prts = cell(1, nv);
end
if isempty(volsel)
    volsel = cell(1, nv);
end
vvol = zeros(1, nv);
for vc = 1:nv
    if ischar(vtcs{vc}) && ...
        exist(vtcs{vc}(:)', 'file') == 2
        loadvtcs(vc) = true;
        try
            vtcs{vc} = xff(vtcs{vc}(:)', 't');
            if ~isxff(vtcs{vc}, 'vtc')
                error( ...
                    'neuroelf:BadArgument', ...
                    'Not a valid VTC filename for run %d.', ...
                    vc ...
                );
            end
        catch ne_eo;
            clearxffobjects(vtcs(loadvtcs));
            clearxffobjects(prts);
            rethrow(ne_eo);
        end
    elseif numel(vtcs{vc}) ~= 1 || ...
       ~isxff(vtcs{vc}, 'vtc')
        clearxffobjects(vtcs(loadvtcs));
        clearxffobjects(prts);
        error( ...
            'neuroelf:BadArgument', ...
            'Cell array must contain either VTC objects or filenames.' ...
        );
    end
    if isempty(volsel{vc})
        volsel{vc} = (1:size(vtcs{vc}.VTCData, 1))';
    else
        volsel{vc} = intersect((1:size(vtcs{vc}.VTCData, 1))', round(volsel{vc}(:)));
    end
    volsel{vc}(isinf(volsel{vc}) | isnan(volsel{vc})) = [];
    vvol(vc) = numel(volsel{vc});
    if ~isempty(xconfound{vc})
        if size(xconfound{vc}, 1) ~= size(vtcs{vc}.VTCData, 1) || ...
            size(xconfound{vc}, 2) >= (0.5 * size(vtcs{vc}.VTCData, 1))
            clearxffobjects(vtcs(loadvtcs));
            clearxffobjects(prts);
            error( ...
                'neuroelf:BadArgument', ...
                'Invalid xconfound size for run %d.', ...
                vc ...
            );
        end
    end
    if isempty(prts{vc})
        if iscell(vtcs{vc}.NameOfLinkedPRT) && ...
           ~isempty(vtcs{vc}.NameOfLinkedPRT)
            prts{vc} = vtcs{vc}.NameOfLinkedPRT{1}(:)';
        elseif ischar(vtcs{vc}.NameOfLinkedPRT)
            prts{vc} = vtcs{vc}.NameOfLinkedPRT(:)';
        end
    end
end

% check headers
vh = vtcs{1};
vs = size(vh.VTCData);
nvol = sum(vvol);
vtr = vh.TR;
for vc = 2:nv
    if vtcs{vc}.Resolution ~= vh.Resolution || ...
        vtcs{vc}.TR ~= vtr || ...
        any(vtcs{vc}.BoundingBox.BBox(:) ~= vh.BoundingBox.BBox(:))
        clearxffobjects(vtcs(loadvtcs));
        clearxffobjects(prts);
        error( ...
            'neuroelf:BadArgument', ...
            'VTCs must match in spatial size/position and temporal settings.' ...
        );
    end
end
vvol = 1 + [0,cumsum(vvol)];

% combine PRTs first
if ~isempty(prts{1})
    for vc = 1:nv
        if ~isempty(prts{vc})
            if lower(prts{vc}.ResolutionOfTime(1)) ~= 'm'
                prts{vc}.ConvertToMS(vtr);
            end
            if volsel{vc}(1) ~= 1
                prts{vc}.ShiftOnsets((1 - volsel{vc}(1)) * vtr);
            end
            if any(diff(volsel{vc}) > 1)
                vss = volsel{vc} - (volsel{vc}(1) - 1);
                vsd = diff(vss);
                vdd = find(vsd > 1);
                for vdc = numel(vdd):-1:1
                    vdt = vss(vdd(vdc)) * vtr;
                    for vcc = 1:numel(prts{vc}.Cond)
                        oo = prts{vc}.Cond(vcc).OnOffsets;
                        if any(oo(:) >= vdt)
                            oo(oo >= vdt) = max(vdt, oo(oo >= vdt) - (vsd(vdd(vdc)) - 1) * vtr);
                            prts{vc}.Cond(vcc).OnOffsets = oo;
                        end
                    end
                end
            end
            for vcc = 1:numel(prts{vc}.Cond)
                oo = prts{vc}.Cond(vcc).OnOffsets;
                vdc = (oo(:, 2) <= 0 | oo(:, 1) == oo(:, 2));
                prts{vc}.Cond(vcc).OnOffsets(vdc, :) = [];
                prts{vc}.Cond(vcc).Weights(vdc, :) = [];
            end
        end
    end
    for vc = 2:nv
        if ~isempty(prts{vc})
            prts{1}.Concatenate(prts{vc}, vvol(vc) - 1, vtr);
        end
    end
    try
        prts{1}.SaveAs(newprt);
    catch ne_eo;
        warning( ...
            'neuroelf:xffError', ...
            'Error saving new PRT: %s.', ...
            ne_eo.message ...
        );
    end
end
clearxffobjects(prts);

% copy first VTC object, and set new NrOfVolumes, VTCData
newvtc = vtcs{1}.CopyObject;
newvtc.NrOfVolumes = nvol;
newvtc.DataType = datatype;
if datatype == 1
    stype = 'uint16';
    mo = uint16(0);
else
    stype = 'single';
    mo = single(0);
    newvtc.FileVersion = 3;
end

% construct transioobject
tobjname = tempname;
tobj = transio(tobjname, 'ieee-le', stype, 0, [nvol, vs(2:4)], 1);
newvtc.VTCData = tobj;

% create temporary object in mem
mo(nvol, vs(2), vs(3)) = 0;

% loop over slices
for sc = 1:vs(4)

    % loop over objects
    for vc = 1:nv

        % get data
        mop = double(vtcs{vc}.VTCData(:, :, :, sc));

        % filter
        if ~isinf(filtsec)
            filtstr = [];
            switch (filttype(1))
                case {'d'}
                    dctnum = 1000 * filtsec / vtcs{vc}.TR;
                    if dctnum <= vtcs{vc}.NrOfVolumes
                        filtstr = struct('tempdct', dctnum);
                    end
                case {'f'}
                    fournum = floor(0.001 * vtcs{vc}.TR * vtcs{vc}.NrOfVolumes / filtsec);
                    if fournum > 0
                        filtstr = struct('tempsc', fournum);
                    end
                case {'p'}
                    polynum = min(15, floor(0.001 * vtcs{vc}.TR * vtcs{vc}.NrOfVolumes / filtsec));
                    if polynum <= vtcs{vc}.NrOfVolumes
                        filtstr = struct('temppoly', polynum);
                    end
            end
            if ~isempty(filtstr)
                if ~isempty(xconfound{vc})
                    filtstr.nuisreg = xconfound{vc};
                end
                mop = tempfilter(mop, filtstr);
            end
        end

        % which transformation
        if transopt == 2
            mop = ztrans(mop, 1);
        elseif transopt == 1
            mop = psctrans(mop, 1);
        end

        % select volumes
        mo(vvol(vc):(vvol(vc + 1)-1), :, :) = mop(volsel{vc}, :, :);
    end

    % store in VTC (transio)
    newvtc.VTCData(:, :, :, sc) = mo;
end

% clear objects
clearxffobjects(vtcs(loadvtcs));

% save VTC under new name
try
    if ~isempty(newprt)
        newvtc.NameOfLinkedPRT = newprt;
    end
    newvtc.SaveAs(targetfile);
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    warning( ...
        'neuroelf:BadArgument', ...
        'Bad target filename given. Please remove ''%s'' manually.', ...
        tobjname ...
    );
    newvtc.ClearObject;
    return;
end

% delete old temp object
delete(tobjname);

% if no argout, clear object
if nargout < 1
    newvtc.ClearObject;
else
    varargout = cell(1, nargout);
    varargout{1} = newvtc;
end
