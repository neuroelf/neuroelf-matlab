function par = readpar(filename, nodata)
% readpar  - read a Philips PAR file
%
% FORMAT:       par = readpar(filename [, nodata])
%
% Input fields:
%
%       filename    PAR filename
%       nodata      flag, if given, do not load transio objects for access
%
% Output fields:
%
%       par         1x1 struct with PAR contents
%
% See also transio

% Version:  v0.9d
% Build:    14061709
% Date:     Jun-17 2014, 9:49 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2014, Jochen Weber
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

% enough arguments ?
if nargin < 1 || ...
   ~ischar(filename) || ...
    isempty(filename) || ...
    exist(filename(:)', 'file') ~= 2
    error( ...
        'neuroelf:BadArguments',...
        'Bad argument or file not found.' ...
    );
end
filename = filename(:)';
if nargin < 2 || ...
   ~islogical(nodata) || ...
    numel(nodata) ~= 1
    nodata = false;
end

% try read
try
    parcont = asciiread(filename);
catch ne_eo;
    error( ...
        'neuroelf:FileNotReadable', ...
        'File not readable: ''%s'' (%s).', ...
        filename, ne_eo.message ...
    );
end

% preliminary parse content
parlines = splittocell(parcont, char([10, 13]), 1, 1);
if numel(parlines) < 30
    error( ...
        'neuroelf:BadFileContent', ...
        'Too few lines in file.' ...
    );
end
for slline = 1:numel(parlines)
    if ~isempty(regexpi(parlines{slline}, 'sl\s+ec\s+dyn'))
        break;
    end
end
if slline == numel(parlines) || ...
    slline < 30
    error( ...
        'neuroelf:BadFileContent', ...
        'Line with ''#sl ec dyn ...'' tokens not found or too early.' ...
    );
end
setlines = parlines(1:(slline - 1));
parlines = parlines((slline + 1):end);
parlines(cellfun('isempty', regexp(parlines, '^\s*\d+'))) = [];

% prepare output
par = struct;
par.PARFile = filename;

% remove comment lines from first part
parset = parsesetlines(grep(setlines, '-vx', '^\#'), '^\.\s+');
par.Parameters = parset;
parnames = fieldnames(parset);

% find settings in par set
noslices = findfield(parnames, {'slice', 'number'});
nodynamics = findfield(parnames, {'dynamic', 'number'});
noechoes = findfield(parnames, {'echo', 'number'});

% check settings
if isempty(noslices) || ...
    isempty(nodynamics)
    error( ...
        'neuroelf:BadFileContent', ...
        'Either number of slices or number of dynamics not found.' ...
    );
end
noslices = parset.(noslices);
nodynamics = parset.(nodynamics);
if ~isempty(noechoes)
    noechoes = parset.(noechoes);
else
    noechoes = 1;
end
nodynsliech = prod([nodynamics, noslices, noechoes]);
if numel(parlines) < nodynsliech
    warning( ...
        'neuroelf:BadFileContent', ...
        'Number of slice description lines mismatch, scan incomplete.' ...
    );
elseif numel(parlines) > nodynsliech
    warning( ...
        'neuroelf:BadFileContent', ...
        'Number of slice description lines mismatch, scan longer than expected.' ...
    );
end

% image definition
imglines = grep(setlines, '-ix', '^#\s+.*(float|integer|string)\)\s*$');
imgset = parsesetlines(imglines, '^#\s+');
par.ImageParameters = imgset;
imgnames = fieldnames(imgset);
imgnum = numel(imgnames);

% build correct representation
imghead = {};
for cc = imgnum:-1:1
    limgname = lower(imgnames{cc});
    limgset = lower(imgset.(imgnames{cc}));
    [li_t{1:3}] = regexp(limgset, '^\((\d+)\*[a-z]+\)$');
    li_t = li_t{3};
    if ~isempty(li_t) && ...
        numel(li_t{1}) == 2
        try
            addnumber = str2double(limgset(li_t{1}(1):li_t{1}(2)));
            imgnum = imgnum + addnumber - 1;
            for ac = addnumber:-1:2
                imghead = [{sprintf('%s_%d', limgname, ac)}; imghead(:)];
            end
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
        end
        imghead = [{[limgname '_1']}; imghead(:)];
    else
        imghead = [{limgname}; imghead(:)];
    end
end
imgnames = imghead;
imgnum = numel(imgnames);

% find settings in img set
sliceno = findfield(imgnames, {'slice', 'number'});
dynamicno = findfield(imgnames, {'dynamic', 'number'});
pixres = findfield(imgnames, {'pixel', 'size'});
imgres = findfield(imgnames, {'resolution'});
imgidx = findfield(imgnames, {'index', 'rec'});

% check settings
if isempty(sliceno) || ...
    isempty(dynamicno) || ...
    isempty(pixres) || ...
    isempty(imgres)
    error( ...
        'neuroelf:BadFileContent', ...
        'Required image header fields not found.' ...
    );
end
sliceno = find(strcmp(imgnames, sliceno));
dynamicno = find(strcmp(imgnames, dynamicno));
pixres = find(strcmp(imgnames, pixres));
imgres = find(strcmp(imgnames, imgres));

% try conversion
valmatrix = u8str2double(gluetostringc(parlines, ';', true));
if size(valmatrix, 1) ~= numel(parlines) || ...
    size(valmatrix, 2) ~= imgnum
    warning( ...
        'neuroelf:BadFileContent', ...
        'Error converting value lines in PAR file.' ...
    );
end
par.MatrixHeaders = imghead;
par.MatrixValues = valmatrix;

% try to locate REC file
[fpath, fname] = fileparts(filename);
if numel(fpath) > 0 && ...
    fpath(end) ~= filesep
    fpath = [fpath filesep];
end
fpath = [fpath fname];
if exist([fpath '.REC'], 'file') == 2
    rfile = [fpath '.REC'];
elseif exist([fpath '.rec'], 'file') == 2
    rfile = [fpath '.rec'];
else
    return;
end
par.RECFile = rfile;

% check pix/img res
if any(diff(valmatrix(:, pixres(1)))) || ...
    any(diff(valmatrix(:, imgres(1)))) || ...
    any(diff(valmatrix(:, imgres(1) + 1)))
    warning( ...
        'neuroelf:FileContNoSupported', ...
        'Unequal image resolutions over time!' ...
    );
    return;
end

% try to guess indexing order
useidx = false;
if ~isempty(imgidx)
    imgidx = find(strcmp(imgnames, imgidx));
    if ~isempty(imgidx)
        [uimgidx, iuimgidx] = unique(valmatrix(:, imgidx(1)));
        if numel(uimgidx) == size(valmatrix, 1)
            useidx = true;
        end
    end
end
par.UseIndex = useidx;

% type of data indexing
if useidx
    dynslc = valmatrix(iuimgidx(:)', [dynamicno(1), sliceno(1)]);
    par.RECDataType = 'Dyn.Slice';
    if diff(valmatrix(iuimgidx(1:2), dynamicno(1))) ~= 0 && ...
        diff(valmatrix(iuimgidx(1:2), sliceno(1))) == 0
        par.RECDataOrder = 'DynFirst';
    elseif diff(valmatrix(iuimgidx(1:2), dynamicno(1))) == 0 && ...
        diff(valmatrix(iuimgidx(1:2), sliceno(1))) ~= 0
        par.RECDataOrder = 'SliceFirst';
    else
        par.RECDataOrder = 'Unknown';
    end
else
    dynslc = valmatrix(:, [dynamicno(1), sliceno(1)]);
    par.RECDataType = 'Slice';
end

% calculate image size and class
pixres = valmatrix(1, pixres(1));
switch (pixres)
    case {8},  dcls = 'uint8';
    case {16}, dcls = 'int16';
    case {32}, dcls = 'single';
    case {64}, dcls = 'double';
    otherwise
        warning( ...
            'neuroelf:FileContNoSupported', ...
            'Pixel bitsize of %d not supported.', ...
            valmatrix(1, pixres(1)) ...
        );
end
pixres = pixres / 8;
imgdims = valmatrix(1, imgres(1):imgres(1)+1);
imgsize = prod(imgdims) * pixres;
par.BytePerPixel = pixres;
par.SliceResolution = imgdims;

% return early
if nodata
    return;
end

% try to create object
try

    % create transio object
    tioobj = struct(transio(rfile, 'ieee-le', dcls, 0, imgdims));
    if useidx
        if size(dynslc, 2) > 1
            for ic = 1:size(dynslc, 1)
                tioobj.IOOffset = (ic - 1) * imgsize;
                par.RECData.Dyn(dynslc(ic, 1)).Slice(dynslc(ic, 2)).IO = ...
                    transio(0, 'makeobject', tioobj);
            end
        else
            for ic = 1:numel(dynslc)
                tioobj.IOOffset = dynslc(ic) * imgsize;
                par.RECData.Slice(ic).IO = ...
                    transio(0, 'makeobject', tioobj);
            end
        end
    else
        if numel(dynslc) > 1
            for ic = 1:numel(dynslc)
                tioobj.IOOffset = (ic - 1) * imgsize;
                par.RECData.Slice(ic).IO = ...
                    transio(0, 'makeobject', tioobj);
            end
        else
            par.RECData.Chunk.IO = ...
                transio(rfile, 'ieee-le', dcls, 0, [dynslc, imgdims]);
        end
    end
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    if isfield(par, 'RECData')
        par = rmfield(par, 'RECData');
    end
end

% %%% subfunctions


function setstr = parsesetlines(setlines, lineprefix)
    setstr = struct;
    for lc = 1:numel(setlines)
        setline = setlines{lc};
        if ~isempty(setline) && ...
            ~isempty(regexpi(setline, lineprefix))
            setline = mstrrep(setline, {lineprefix, '\s+'}, {'', ' '}, 1);
            colonpos = find(setline == ':');
            if numel(colonpos) == 1 && ...
                colonpos > 1
                labline = makelabel(deblank(regexprep(setline(1:colonpos - 1), ...
                    '[^a-zA-Z\s\.].*$', '')));
                valline = mstrrep(setline(colonpos + 1:end), ...
                    {'^\s+', '\s+$'}, {'', ''}, 1);
            else
                setlinec = splittocell(setline, char([9, 32]), 1, 1);
                if numel(setlinec) < 2
                    continue;
                end
                labline = makelabel(deblank(regexprep(gluetostring( ...
                    setlinec(1:end-1), ' '), '[^a-zA-Z_\s\.].*$', '')));
                valline = setlinec{end};
            end
            if numel(labline) < 2 || ...
                (any('vV' == labline(1)) && labline(2) == '_')
                continue;
            end
            vallineok = regexp(valline, '^[0-9\+\-e\.\s]+$');
            if ~isempty(vallineok) && ...
               (isa(vallineok, 'double') || ...
                (iscell(vallineok) && ~isempty(vallineok{1})))
                try
                    valline = eval(['[' valline ']']);
                catch ne_eo;
                    neuroelf_lasterr(ne_eo);
                    valline = 0;
                end
            end
            setstr.(labline) = valline;
        end
    end
% end of function setstr = parsesetlines(setlines)

function fname = findfield(fnames, parts)
    lfnames = lower(fnames);
    if ~iscell(parts)
        parts = {parts};
    end
    parts = lower(parts);
    pnum = numel(parts);
    for fc = 1:numel(fnames)
        fmatch = true;
        for pc = 1:pnum
            if isempty(strfind(lfnames{fc}, parts{pc}))
                fmatch = false;
                break;
            end
        end
        if fmatch
            fname = fnames{fc};
            return;
        end
    end
    fname = '';
% end of function fname = findfield(fnames, parts)
