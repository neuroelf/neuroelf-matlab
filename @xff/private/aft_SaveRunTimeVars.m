function xo = aft_SaveRunTimeVars(xo, savecomp)
% AnyFileType::SaveRunTimeVars  - updates the RunTimeVars field on disk
%
% FORMAT:       object.SaveRunTimeVars;
%
% No input / output fields.
%
% TYPES: ALL
%
% Using: catstruct, emptystruct, mainver, renamefile.

% Version:  v1.1
% Build:    16022911
% Date:     Feb-29 2016, 11:23 AM EST
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

% global neuroelf methods
global ne_methods;

% only valid for single file
if numel(xo) ~= 1 || ~xffisobject(xo, true) || xo.L(1) == 'X'
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
if nargin < 2 || ~islogical(savecomp) || numel(savecomp) ~= 1
    savecomp = true;
end

% check filename
if isempty(xo.F)
    error('neuroelf:xff:badFilename', 'File not yet saved. Use SaveAs method first.');
end

% don't allow volume marker
if ~isempty(regexpi(xo.F, ',\d+$'))
    error('neuroelf:xff:badArgument', 'Saving of sub-volumes not permitted.');
end

% see if run time vars already exists
bc = xo.C;
RunTimeVars = struct;

% support for gziped files
if isfield(xo.H, 'GZIPext') && ischar(xo.H.GZIPext) && strcmpi(xo.H.GZIPext, '.gz') && ...
    isfield(xo.H, 'GZIPfile') && ischar(xo.H.GZIPfile) && ~isempty(xo.H.GZIPfile)
    [filenp, filenn, filene] = fileparts(xo.H.GZIPfile);
else
    [filenp, filenn, filene] = fileparts(xo.F);
end
mmversion = ne_methods.mainver();
try
    filenm = fopen([filenp '/' filenn '.rtv']);
    if filenm > 0
        fclose(filenm);
        if mmversion >= 5
            filenm = load([filenp '/' filenn '.rtv'], '-mat');
        else
            filenm = load('-mat', [filenp '/' filenn '.rtv']);
        end
        if ~isfield(filenm, 'RunTimeVars') || ~isstruct(filenm.RunTimeVars) || ...
            numel(filenm.RunTimeVars) ~= 1
            error('neuroelf:xff:invalidFile', ...
                'The RunTimeVars .mat file cannot be updated.');
        end
        RunTimeVars = filenm.RunTimeVars;
    end
catch xfferror
    neuroelf_lasterr(xfferror);
end

% update RunTimeVars
RunTimeVars.(filene(2:end)) = bc.RunTimeVars;

% gather additional content
bcf = fieldnames(xo.C);
for bcfc = 1:numel(bcf)
    if isstruct(bc.(bcf{bcfc})) && isfield(bc.(bcf{bcfc}), 'RunTimeVars') && ...
       ~isempty(bc.(bcf{bcfc})) && isstruct(bc.(bcf{bcfc})(1).RunTimeVars) && ...
        numel(bc.(bcf{bcfc})(1).RunTimeVars) == 1
        rtvsc = {bc.(bcf{bcfc}).RunTimeVars};
        if numel(rtvsc) > 1
            rtvsc = rtvsc(:);
            rtvsf = fieldnames(rtvsc{1});
            for rtvscc = 2:numel(rtvsc)
                if ~isstruct(rtvsc{rtvscc}) || numel(rtvsc{rtvscc}) ~= 1
                    rtvsc{rtvscc} = ne_methods.emptystruct(rtvsf, [1, 1]);
                end
            end
            RunTimeVars.(filene(2:end)).(bcf{bcfc}) = ne_methods.catstruct(rtvsc{:});
        else
            RunTimeVars.(filene(2:end)).(bcf{bcfc}) = rtvsc{1};
        end
    end
end

% try to save to file
try
    rtvfname = sprintf('%08x', floor(2^32 * rand(1, 4)));
    if mmversion < 5
        if savecomp
            save('-v7', [filenp '/' filenn '_' rtvfname '_rtv.mat'], 'RunTimeVars');
        else
            save('-v6', [filenp '/' filenn '_' rtvfname '_rtv.mat'], 'RunTimeVars');
        end
    elseif mmversion < 7 || savecomp
        save([filenp '/' filenn '_' rtvfname '_rtv.mat'], 'RunTimeVars');
    else
        save([filenp '/' filenn '_' rtvfname '_rtv.mat'], 'RunTimeVars', '-v6');
    end
    try
        ne_methods.renamefile([filenp '/' filenn '_' rtvfname '_rtv.mat'], ...
            [filenp '/' filenn '.rtv']);
    catch xfferror
        neuroelf_lasterr(xfferror);
        warning('neuroelf:xff:fileWriteError', ...
            'Error renaming temporary RunTimeVars file to final location.');
    end

% and report errors
catch xfferror
	rethrow(xfferror);
end

% if all went well, set handle
xo.H.RunTimeVarsSaved = true;
