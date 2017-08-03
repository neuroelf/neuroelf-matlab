function [varargout] = renamefile(varargin)
% renamefile  - renames a file to another filename
%
% FORMAT:       renamefile(FromFile, ToFileOrDir);
%
% Input fields:
%
%       FromFile    filename(s) of files being renamed
%                   either a char array for single file or
%                   cell array with multiple char arrays
%       ToFileOrDir one of the following:
%                   - single filename as target filename
%                   - single foldername as target for one or
%                     multiple files
%                   - cell array with multiple targetnames
%                     (dims must match with input array)
%
% Note: there is a special usage you can use to replace a pattern
%       in one or multiple filenames:
%
% renamefile(FileOrListOfFiles, FromPattern, ToPattern);

% Version:  v0.9c
% Build:    11071411
% Date:     Jul-14 2011, 11:26 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, Jochen Weber
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
if nargin < 2
    error( ...
        'neuroelf:BadArgument', ...
        'Too few arguments. Try ''help %s''.', ...
        mfilename ...
    );
end

% act on special case
if nargin > 2 && ...
    ischar(varargin{2}) && ...
   ~isempty(varargin{2}) && ...
    ischar(varargin{3})
    if ischar(varargin{1}) || ...
       (iscell(varargin{1}) && ...
       ~isempty(varargin{1}) && ...
        ischar(varargin{1}{1}))
        try
            rval = renamefile(varargin{1}, ...
                       strrep(varargin{1}, varargin{2}, varargin{3}));
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            rval = -1;
        end
    else
        error( ...
            'neuroelf:BadArgument', ...
            'Bad argument given.' ...
        );
    end
    if nargout > 0
        varargout{1} = rval;
    end
    return;
end

FromFile = varargin{1};
ToFileOrFolder = varargin{2};
if nargout > 0
    varargout{1}=-1;
end
if iscell(FromFile) && ...
    numel(FromFile) == 1
    FromFile = FromFile{1};
end
if iscell(ToFileOrFolder) && ...
    numel(ToFileOrFolder) == 1
    ToFileOrFolder = ToFileOrFolder{1};
end

if (~ischar(FromFile) && ...
    ~iscell(FromFile)) || ...
   (~ischar(ToFileOrFolder) && ...
    ~iscell(ToFileOrFolder))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad argument.' ...
    );
end
if ischar(FromFile) && ...
   ~ischar(ToFileOrFolder)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad argument.' ...
    );
end

if ischar(FromFile)
    if exist(FromFile, 'file') < 2
        if nargout > 0
            varargout{1} = 8;
        end
        warning( ...
            'neuroelf:RenameFailed', ...
            'Source file ''%s'' not found.', ...
            FromFile ...
        );
        return;
    end
    if exist(ToFileOrFolder, 'dir') == 7
        TargetFile = i_movetodir(FromFile, ToFileOrFolder);
        if exist(TargetFile, 'file') < 2
            if nargout > 0
                varargout{1} = 16;
            end
            warning( ...
                'neuroelf:RenameFailed', ...
                'Couldn''t move source to destination folder.' ...
            );
            return;
        end
    else
        TargetFile = i_movetofile(FromFile, ToFileOrFolder);
        if exist(TargetFile, 'file') < 2
            if nargout > 0
                varargout{1} = 16;
            end
            warning( ...
                'neuroelf:RenameFailed', ...
                'Couldn''t move source to destination file.' ...
            );
            return;
        end
    end

else

    % initialize counters
    NumOfFiles = numel(FromFile);
    CompleteSuccess = 1;

    if (ischar(ToFileOrFolder) && ...
        exist(ToFileOrFolder, 'dir') ~= 7) || ...
        (iscell(ToFileOrFolder) && ...
         numel(ToFileOrFolder) ~= NumOfFiles)
        error( ...
            'neuroelf:BadArgument', ...
            'Bad sized argument.' ...
        );
    end

    if iscell(ToFileOrFolder)
        for FileNumber = 1:NumOfFiles
            if exist(ToFileOrFolder{FileNumber}, 'dir') == 7
                TargetFile = ...
                    i_movetodir(FromFile{FileNumber}, ToFileOrFolder{FileNumber});
            else
                TargetFile = ...
                    i_movetofile(FromFile{FileNumber}, ToFileOrFolder{FileNumber});
            end
            if exist(TargetFile, 'file') < 2
                CompleteSuccess = 0;
                warning( ...
                    'neuroelf:RenameFailed', ...
                    'Couldn''t rename ''%s'' to ''%s''.', ...
                    FromFile{FileNumber}, TargetFile ...
                );
            end
        end
    else
        for FileNumber = 1:NumOfFiles
            TargetFile = i_movetodir(FromFile{FileNumber}, ToFileOrFolder);
            if exist(TargetFile, 'file') < 2
                CompleteSuccess = 0;
                warning( ...
                    'neuroelf:RenameFailed', ...
                    'Couldn''t rename ''%s'' to ''%s''.', ...
                    FromFile{FileNumber}, TargetFile ...
                );
            end
        end
    end

    if nargout > 0
        if CompleteSuccess ~= 1
            varargout{1} = 1;
        else
            varargout{1} = 0;
        end
    end

end

if nargout > 0 && ...
    varargout{1} == -1
    varargout{1} = 0;
end



% sub functions



function tgf = i_movetodir(ff, tf)
    if any(ff == filesep)
        [sp{1:3}] = fileparts(ff);
        sn = [sp{2} sp{3}];
    else
        sn = ff;
    end
    tgf = [tf filesep sn];
    if ispc
        if exist('movefile', 'builtin') == 5
            movefile(ff, tgf);
        else
            system(['move "' ff '" "' tgf '"']);
        end
    else
        system(['mv   "' ff '" "' tgf '"']);
    end
% end of function tgf = i_movetodir(ff, tf)

function tgf = i_movetofile(ff, tf)
    if ~any(tf == filesep) && ...
        any(ff == filesep)
        sp = fileparts(ff);
        tgf = [sp filesep tf];
    else
        tgf = tf;
    end
    if ispc
        ff = strrep(ff, '/', '\');
        tgf = strrep(tgf, '/', '\');
    end
    if exist('movefile', 'builtin') == 5
        movefile(ff, tf);
    else
        if ispc
            system(['move "' ff '" "' tgf '"']);
        else
            system(['mv   "' ff '" "' tgf '"']);
        end
    end
% end of function tgf = i_movetofile(ff, tf)
