function [folder, file, ext] = mfileparts(filenames)
% mfileparts  - run fileparts along multiple files
%
% FORMAT:       [folders, files, exts] = mfileparts(filenames)
%
% Input fields:
%
%       filenames   cell array with filenames
%
% Output fields:
%
%       folders     cell array with folder name(s)
%       files       cell array with filenames
%       exts        cell array with extensions

% Version:  v0.9c
% Build:    13110615
% Date:     Nov-06 2013, 3:55 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2013, Jochen Weber
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
   (~ischar(filenames) && ...
    ~iscell(filenames))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if ischar(filenames)
    filenames = cellstr(filenames);
else
    filenames = filenames(:);
    if ~all(cellfun(@ischar, filenames))
        error( ...
            'neuroelf:BadArgument', ...
            'Bad or missing argument.' ...
        );
    end
end

% create main output
folder = cell(numel(filenames), 1);

% not all outputs required
if nargout < 3

    % only one output required
    if nargout < 2

        % iterate
        for fc = 1:numel(folder)
            folder{fc} = fileparts(filenames{fc});
        end

        % all folders the same?
        if all(strcmp(folder, folder{1}))
            folder = folder(1);
        end

    % two outputs required
    else

        % create second output
        file = cell(numel(folder), 1);

        % iterate
        for fc = 1:numel(folder)
            [folder{fc}, file{fc}] = fileparts(filenames{fc});
        end

        % all folders the same?
        if all(strcmp(folder, folder{1}))
            folder = folder(1);
        end
    end

% three outputs
else

    % create second and third output
    file = cell(numel(folder), 1);
    ext = cell(numel(file), 1);

    % iterate
    for fc = 1:numel(folder)
        [folder{fc}, file{fc}, ext{fc}] = fileparts(filenames{fc});
    end

    % all folders/extensions the same?
    if all(strcmp(folder, folder{1}))
        folder = folder(1);
    end
    if all(strcmp(ext, ext{1}))
        ext = ext(1);
    end
end
