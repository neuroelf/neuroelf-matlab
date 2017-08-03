function result = xffcompare(newfile, compfile)
% xffcompare  - compares two xff objects (or folders with files)
%
% FORMAT:       result = xffcompare(newfile, compfile)
%
% Input fields:
%
%       newfile     file (or folder) that has to be compared
%       compfile    reference file (or folder with) copy (copies)
%
% Output fields:
%
%       result      text containing results of comparison
%
% Note: possible results are
%       FNFND - new file not found
%       RMISS - reference file is missing
%       ERROR - error reading new or reference file
%       IDENT - files are identical
%       TDIFF - xff types are different
%       FDIFF - types are the same, but fields are different
%       VDIFF - simple version differences (loading outcome is equal)
%       CDIFF - content is different (different fields will be listed)
%       for different fields, the differences can be
%       CLASD - classes differ
%       DIMSD - dimensions differ
%       VALUE - for 1x1 fields, value is not equal
%       MINOR - minor differences (e.g. for text fields)
%       SUBTL - subtle differences (e.g. rounding errors)
%       MAJOR - major content wise differences
%       in case of major differences, a quantification is given in
%       EQELM (equal elements ratio), MABSD (mean of absolute deviation),
%       STDEV (standard deviation)

% Version:  v1.1
% Build:    16012314
% Date:     Jan-23 2016, 2:00 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2016, Jochen Weber
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
   ~ischar(newfile) || ...
    isempty(newfile) || ...
   ~ischar(compfile) || ...
    isempty(compfile)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
newfile = newfile(:)';
compfile = compfile(:)';

% file or folder
nft = exist(newfile, 'file');
cft = exist(compfile, 'file');
isf = 0;
if nft == 2
    if cft == 2
        isf = 2;
    elseif cft == 7
        [nfp] = fileparts(newfile);
        compfile = [compfile '/' nfp{2} nfp{3}];
        if exist(compfile, 'file') == 2
            isf = 2;
        end
    end
elseif nft == 7
    if cft == 2
        [cfp] = fileparts(compfile);
        newfile = [newfile '/' cfp{2} cfp{3}];
        if exist(newfile, 'file') == 2
            isf = 2;
        end
    elseif cft == 7
        isf = 7;
    end
else
    result = 'FNFND';
    return;
end
if isf == 0
    result = 'RMISS';
    return;
end

% for folders, iteratively call for all source file
if isf == 7
    nfiles = findfiles(newfile, '*.*', 'relative=');
    rfiles = cell(1, numel(nfiles));
    for nfc = 1:numel(nfiles)
        cfile = [compfile '/' nfiles{nfc}];
        if exist(cfile, 'file') ~= 2
            rfiles{nfc} = 'RMISS';
        else
            rfiles{nfc} = xffcompare([newfile '/' nfiles{nfc}], cfile);
        end
        rfiles{nfc} = sprintf('%s: %s\n', nfiles{nfc}, rfiles{nfc});
    end
    result = gluetostring(rfiles, '');
    return;
end

% compare file sizes
try
    nfsize = filesize(newfile);
    cfsize = filesize(compfile);
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    result = 'ERROR';
    return;
end

% if identical in size, simply try binary match
if nfsize == cfsize

    % for small files, use entire file with binread
    if nfsize <= 16777216
        try
            if all(binread(newfile) == binread(compfile))
                result = 'IDENT';
                return;
            end
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            return;
        end

    % for larger files
    else

        % open both files
        try
            nfid = 0;
            cfid = 0;
            nfid = fopen(newfile);
            if nfid < 1
                error('ERROR_FOPEN');
            end
            cfid = fopen(compfile);
            if cfid < 1
                error('ERROR_FOPEN');
            end

            % and read in 16MB chunks
            while nfsize > 16777216
                if any(fread(nfid, [1, 16777216], '*uint8') ~= ...
                        fread(cfid, [1, 16777216], '*uint8'))
                    break;
                end
                nfsize = nfsize - 16777216;
            end
            if nfsize <= 16777216
                if all(fread(nfid, [1, nfsize], '*uint8') == ...
                        fread(cfid, [1, nfsize], '*uint8'))
                    nfsize = 0;
                end
            end

            % close files
            fclose(nfid);
            fclose(cfid);

            % if result reached
            if nfsize == 0
                result = 'IDENT';
                return;
            end

        % on error
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            if nfid > 0
                try
                    fclose(nfid);
                catch ne_eo;
                    neuroelf_lasterr(ne_eo);
                end
            end
            if cfid > 0
                try
                    fclose(cfid);
                catch ne_eo;
                    neuroelf_lasterr(ne_eo);
                end
            end
            result = 'ERROR';
            return;
        end
    end
end

% files are not identical (either in size or binary content)
% try loading both as xff objects with transio access
xffroot = xff();
tiostr = xffroot.TransIOSize(8192);
xffobj = cell(1, 2);
try
    xffobj{1} = xff(newfile);
    xffobj{2} = xff(compfile);
    if ~isxff(xffobj{1}) || ...
       ~isxff(xffobj{2})
        error('READERROR');
    end

% on error, give back ERROR
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    clearxffobjects(xffobj);
    result = 'ERROR';
    xffroot.TransIOSize(tiostr);
    return;
end

% get content and clear objects
ntp = xffobj{1}.Filetype;
ctp = xffobj{2}.Filetype;
nfo = getcont(xffobj{1});
cfo = getcont(xffobj{2});
clearxffobjects(xffobj);
xffroot.TransIOSize(tiostr);

% check for same type
if ~strcmp(ntp, ctp)
    result = 'TDIFF';
    return;
end

% check for same fields
nff = fieldnames(nfo);
cff = fieldnames(cfo);
if numel(nff) ~= numel(cff) || ...
    ~all(strcmp(nff, cff))
    result = 'FDIFF';
    return;
end

% check each field (suppose equality)
result = 'VDIFF';
fdiff = cell(1, numel(nff));
for fc = 1:numel(nff)
    fdiff{fc} = '';
    fcls = lower(class(nfo.(nff{fc})));
    if ~strcmpi(fcls, class(cfo.(nff{fc})))
        fdiff{fc} = [nff{fc} ':CLASD,'];
        continue;
    end
    fsiz = size(nfo.(nff{fc}));
    if ~isequal(fsiz, size(cfo.(nff{fc})))
        fdiff{fc} = [nff{fc} ':DIMSD,'];
        continue;
    end
    if isequal(fsiz, [1, 1])
        if ~isequal(nfo.(nff{fc}), cfo.(nff{fc}))
            fdiff{fc} = [nff{fc} ':VALUE,'];
        end
        continue;
    end
    switch (fcls)
        case {'char'}
            fdiff{fc} = [nff{fc} ':MINOR,'];
        case {'transio'}
        otherwise
            fdiff{fc} = [nff{fc} ':UNDEF,'];
    end
end
fdiff = gluetostring(fdiff, '');
if ~isempty(fdiff)
    result = ['CDIFF (' fdiff(1:end-1) ')'];
end
