function xo = aft_SaveAs(xo, newfile, dtitle)
% AnyFileType::SaveAs  - method for any xff type
%
% FORMAT:       obj.SaveAs([newfilename, dtitle]);
%
% Input fields:
%
%       newfilename save-as filename, if not give, UI-based
%       dtitle      if given, override UI default title
%
% No output fields:
%
%       obj         xff object with newly set filename
%
% TYPES: ALL
%
% Using: bffio, tffio.

% Version:  v1.1
% Build:    16040820
% Date:     Apr-08 2016, 8:01 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2015, 2016, Jochen Weber
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

% global neuroelf methods and object list (for filename update)
global ne_methods xffsngl;

% check arguments
if nargin < 1 || numel(xo) ~= 1 || ~xffisobject(xo, true) || ~isfield(xo.S, 'FFTYPE')
    error('neuroelf:xff:badArguments', 'Invalid call to %s.', mfilename);
end

% get objects super-struct
scf = xo.F;

% arguments
ndtitle = '';
if nargin > 1 && ischar(newfile) && ~any(newfile == '.')
    ndtitle = newfile(:)';
    newfile = '';
end
if nargin < 2 || ~ischar(newfile) || isempty(newfile)
    if isempty(scf)
        newfile = ['*.' xo.S.Extensions{1}];
    else
        [of{1:3}] = fileparts(scf);
        newfile = ['*' of{3}];
    end
    if strcmpi(xo.S.Extensions{1}, 'hdr') && strcmpi(xo.C.FileMagic, 'n+1')
        newfile = regexprep(newfile, '^(.*)\.hdr$', '$1.nii', 'preservecase');
    end
end
if nargin < 3 || ~ischar(dtitle) || isempty(dtitle)
    if ~isempty(ndtitle)
        dtitle = ndtitle;
    else
        dtitle = 1;
    end
end
newfile = newfile(:)';

% don't allow volume marker
if ~isempty(regexpi(newfile, ',\d+$'))
    error('neuroelf:xff:badArgument', 'Saving of sub-volumes not permitted.');
end

% make absolute name
[isabs{1:2}] = isabsolute(newfile);
newfile = isabs{2};
[fnparts{1:3}] = fileparts(newfile);

% check for "*.???"
if isempty(fnparts{2})
    fnparts{2} = '*';
end
if isempty(fnparts{3}) || strcmp(fnparts{3}, '.')
    fnparts{3} = ['.' xo.S.Extensions{1}];
end
if any(fnparts{2} == '*')
    extensions   = xffsngl.EXT;
    file_formats = xffsngl.FF;
    filename = xffrequestfile(dtitle, ...
        [fnparts{1} filesep fnparts{2} fnparts{3}], extensions, file_formats, true);
    if isempty(filename)
        return;
    end
    if iscell(filename)
        filename = filename{1};
    end
    [isabs{1:2}] = isabsolute(filename);
    newfile = isabs{2};
end

% what to do
try
    switch (lower(xo.S.FFTYPE))
        case 'bff'
            xo.C = ne_methods.bffio(newfile, xo.S, xo.C);
            xo.F = newfile;
        case 'tff'
            [xo.C, xo.F] = ne_methods.tffio(newfile, xo.S, xo.C);
        otherwise
            error('neuroelf:xff:invalidFileType', ...
                'Type not recognized (?FF): %s.', xo.S.FFTYPE);
    end
catch xfferror
    error('neuroelf:xff:errorSavingFile', 'Error saving file %s: %s.', xo.F(:)', xfferror.message);
end

% update filename
olup = ne_methods.findfirst(strcmpi(xffsngl.OBJS(:, 3), xo.L));
if ~isempty(olup)
    xffsngl.OBJS{olup, 1} = xo.F;
end

% remove GZIP temp file?
if isfield(xo.H, 'GZIPext') && ischar(xo.H.GZIPext) && strcmpi(xo.H.GZIPext, '.gz') && ...
    isfield(xo.H, 'GZIPfile') && ischar(xo.H.GZIPfile) && ~isempty(xo.H.GZIPfile) && ...
    exist([xo.H.GZIPfile xo.H.GZIPext], 'file') == 2 && ~isempty(scf) && exist(scf, 'file') == 2
    try
        delete(scf);
    catch xfferror
        neuroelf_lasterr(xfferror);
    end
end
        
% remove GZIPext/file from handles
if isfield(xo.H, 'GZIPext')
    xo.H = rmfield(xo.H, 'GZIPext');
end
if isfield(xo.H, 'GZIPfile')
    xo.H = rmfield(xo.H, 'GZIPfile');
end

% then see if RunTimeVars are to be saved as well
if isfield(xo.C.RunTimeVars, 'AutoSave') && islogical(xo.C.RunTimeVars.AutoSave) && ...
    numel(xo.C.RunTimeVars.AutoSave) == 1 && xo.C.RunTimeVars.AutoSave

    % try automatic saving
    try
        aft_SaveRunTimeVars(xo);
    catch xfferror
        warning('neuroelf:xff:errorSavingFile', ...
            'Error saving RunTimeVars file: %s.', xfferror.message);
    end
end
