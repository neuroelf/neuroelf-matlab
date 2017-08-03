function [o_status, o_message] = mkadir(varargin)
% mkadir  - create an absolute directory (unlike to MATLAB's mkdir)
%
% FORMAT:         mkadir <dirname> [-p]
%
%    dirname      string absolute directory to create
%    -p           when given, try to create parent dir(s) first
%
% the usage is straight forward...
%
% mkadir c:\temp\this
% mkadir c:\temp\this\is\just\a\test -p
%
% st = mkadir(['/tmp/' temppath '/.log');
%
% See also mkdir.

% Version:  v0.9a
% Build:    10051716
% Date:     May-17 2010, 10:48 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, Jochen Weber
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
if nargin < 1
    error( ...
        'neuroelf:TooFewArguments',...
        'Too few arguments. Try ''help %s''.',...
        mfilename ...
    );
end

% sanity checks
for acount = 1:nargin
    if ~ischar(varargin{acount})
        error( ...
            'neuroelf:BadArgument',...
            'All arguments must be of type char.' ...
        );
    end
end

% extract parent for MATLAB's mkdir and check if mkadir will work
dirname = varargin{1};
if dirname(end) == filesep
    dirname(end) = [];
end

% get parts
[parent, name, ext] = fileparts(dirname);

% if parent does not yet exist
if exist(parent,'dir') ~= 7

    % only if -p argument is given
    if nargin > 1 && ...
        ischar(varargin{2}) && ...
        numel(varargin{2}) == 2 && ...
        all(lower(varargin{2}(:)') == '-p')

        % call recursively
        status = mkadir(parent,'-p');

        % status check
        if status < 1
            error( ...
                'neuroelf:DirNotCreated',...
                'Couldn''t create dir ''%s''.',...
                strrep(parent,'\','\\') ...
            );
        end

    % otherwise
    else

        % give error
        error( ...
            'neuroelf:FolderNotExists',...
            'Parent folder ''%s'' doesn''t exist. Try ''%s <DIR> -p''.',...
            strrep(parent,'\','\\'),...
            mfilename ...
        );
    end
end

% reassemble wanted filename
if ~isempty(ext)
    name = [name ext];
end

% let MATLAB to the rest
[status, message] = mkdir(parent, name);

% handle return values correctly
if status < 1
    error( ...
        'neuroelf:DirNotCreated',...
        'Couldn''t create ''%s'' in ''%s''.',...
        name,...
        strrep(parent,'\','\\') ...
    );
end

% handle output
if nargout > 0
    o_status  = status;
end
if nargout > 1
    o_message = message;
end
