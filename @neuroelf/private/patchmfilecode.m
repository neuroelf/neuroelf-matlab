function patchmfilecode(sourcefile, targetfile, patches, opts)
%PATCHMFILECODE  Patch an M-file by applying changes and store a new file.
%   PATCHMFILECODE(SOURCEFILE, TARGETFILE, PATCHES) reads the M-file in
%   SOURCEFILE, applies the patches specified in PATCHES and stores the
%   resulting output in TARGETFILE.
%
%   Both SOURCEFILE *and* TARGETFILE must be specified as full path names.
%
%   PATCHES is a Px3 cell array, where the first column is used to instruct
%   the function what kind of patch to apply; valid options are
%
%   'regexprep'     - applies a "regexprep" using the second/third columns
%   'replacelines'  - replaces the lines specified (in the second column)
%                     with the content in the third column
%   'strrep'        - applies a "strrep" using the second/third columns
%
%   PATCHMFILECODE(SOURCEFILE, TARGETFILE, PATCHES, OPTS) allows to also
%   specify the following options by passing a 1x1 struct OPTS with fields
%
%       .overwrite  overwrite existing TARGETFILE (otherwise, error out)
%
%   If TARGETFILE exists, an error will be thrown, unless the .overwrite
%   option is set to (1x1 logical) true.
%
%   Example:
%
%   PATCHMFILECODE([neuroelf_path '/neuroelf_version.m'], ...
%       [neuroelf_path '/neuroelf_version_patched.m'], ...
%       {'strrep', 'Jochen Weber', 'The NeuroElf'}, struct('overwrite', true));

% Version:  v1.1
% Build:    16060911
% Date:     Jun-09 2016, 11:05 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2016, Jochen Weber
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
if nargin < 3 || ~ischar(sourcefile) || ~ischar(targetfile) || ~iscell(patches) || ...
    isempty(sourcefile) || isempty(targetfile) || isempty(patches) || size(patches, 2) ~= 3 || ...
    exist(sourcefile(:)', 'file') ~= 2 || any(cellfun('isempty', patches(:, 1))) || ...
    any(cellfun('isempty', patches(:, 2))) || ~all(cellfun(@ischar, patches(:, 1))) || ...
    any(cellfun('isempty', regexpi(patches(:, 1), '^(regexprep|replacelines|strrep)$')))
    error('neuroelf:general:badArgument', 'Bad or missing argument.');
end
sourcefile = sourcefile(:)';
targetfile = targetfile(:)';
spath = fileparts(sourcefile);
tpath = fileparts(targetfile);
if isempty(spath) || isempty(tpath)
    error('neuroelf:general:badArgument', 'Both sourcefile and targetfile must be absolute paths.');
end
if nargin < 4 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'overwrite') || ~islogical(opts.overwrite) || numel(opts.overwrite) ~= 1
    opts.overwrite = false;
end
if ~opts.overwrite && exist(targetfile, 'file') == 2
    error('neuroelf:badOption:noFileOverwrite', 'Target file exists and overwrite not requested.');
end

% try to read source file
try
    sfid = -1;
    sfid = fopen(sourcefile, 'r');
    if sfid < 1
        error('neuroelf:fileError:fileOpenError', 'Error opening source file.');
    end
    filecont = fread(sfid, [1, Inf], 'char=>char');
catch readerror
    if sfid > 0
        fclose(sfid);
    end
    rethrow(readerror);
end

% apply patches
for pc = 1:size(patches, 1)

    % patch type, from and to specs
    ptype = lower(patches{pc, 1}(3)); 
    fromspec = patches{pc, 2};
    tospec = patches{pc, 3};

    % regexprep
    if ptype == 'g'
        if ~ischar(tospec)
            error('neuroelf:general:badArgument', 'Invalid replace-to spec in patch %d.', pc);
        end

        % apply
        filecont = regexprep(filecont, fromspec(:)', tospec(:)');

    % strrep
    elseif ptype == 'r'
        if ~ischar(tospec)
            error('neuroelf:general:badArgument', 'Invalid replace-to spec in patch %d.', pc);
        end

        % apply
        filecont = strrep(filecont, fromspec(:)', tospec(:)');

    % replace lines
    else

        % find lines breaks
        lb = regexp(filecont, '(\r\n|\r|\n)');
        lb(end+1) = numel(filecont) + 1;
        numlines = numel(lb);

        % detect kind of linebreak
        if numlines < 2
            lbc = char(10);
        elseif numlines == 2
            if lb == numel(filecont)
               lbc = filecont(end); 
            else
                if all(filecont(lb(1):lb(1)+1) == char([13, 10]))
                    lbc = char([13, 10]);
                else
                    lbc = filecont(lb(1));
                end
            end
        elseif all(filecont(lb(1):lb(1)+1) == char([13, 10]))
            lbc = char([13, 10]);
        else
            lbc = filecont(lb(1));
        end

        % needs contiguous line numbers
        if ~isa(fromspec, 'double') || isempty(fromspec) || ...
            any(isinf(fromspec(:)) | isnan(fromspec(:)) | fromspec(:) < 1 | fromspec(:) > numlines) || ...
            any(fromspec(:) ~= round(fromspec(:))) || numel(fromspec) ~= numel(unique(fromspec(:))) || ...
            (max(fromspec(:)) - min(fromspec(:))) ~= (numel(fromspec) - 1)
            error('neuroelf:general:badArgument', 'Invalid replacelines spec in patch %d.', pc);
        end
        fromline = min(fromspec(:));
        toline = max(fromspec(:)) + 1;

        % compile replacement
        if isempty(tospec) || (iscell(tospec) && numel(tospec) == 1 && isempty(tospec{1}))
            tospec = '';
        elseif ischar(tospec)
            tospec = tospec(:)';
            if ~any(tospec(end) == char([10, 13]))
                tospec = [tospec, lbc];
            end
            tospec = regexprep(tospec, '(\r\n|\r|\n)', lbc);
        elseif iscell(tospec) && all(cellfun(@ischar, tospec(:)))
            tospec = [tospec(:)'; repmat({lbc}, 1, numel(tospec))];
            tospec = sprintf('%s', tospec{:});
        else
            error('neuroelf:general:badArgument', 'Invalid replace-to specification in patch %d.', pc);
        end

        % replace
        filecont = [filecont(1:lb(fromline)-1), tospec, filecont(lb(toline):end)];
    end
end

% write out altered file content
try
    tfid = -1;
    tfid = fopen(targetfile, 'w');
    if tfid < 1
        error('neuroelf:fileError:fileOpenError', 'Error opening target file.');
    end
    fwrite(tfid, filecont, 'char');
    fclose(tfid);
catch writeerror
    if tfid > 0
        fclose(tfid);
    end
    rethrow(writeerror);
end
