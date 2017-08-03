function [isabs, abspath] = isabsolute(tpath)
% isabsolute  - returns true if the given path is absolute
%
% FORMAT:       isabs = isabsolute(testpath)
%
% Input fields:
%
%       testpath    string to test for absolute path content
%
% Output fields:
%
%       isabs       true for absolute paths, false for relative
%
% See also EXIST
%
% Using: ostype.

% Version:  v0.9d
% Build:    14082216
% Date:     Aug-22 2014, 4:45 PM EST
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

% persistent ostype
persistent machine;
if isempty(machine) || ...
   ~ischar(machine)
    try
        nelf = neuroelf;
        ostype = nelf.ostype();
        machine = ostype.machine;
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        if ispc
            machine = 'WIN';
        else
            machine = 'UNIX';
        end
    end
end

% argument check
if nargin < 1 || ...
   ~ischar(tpath) || ...
    isempty(tpath)
    error( ...
        'neuroelf:BadArgument', ...
        'One non-empty 1xN char path argument is required.' ...
    );
end

% make valid char argument
tpath = tpath(:)';

% default is false
isabs = false;
abspath = tpath;

% machine dependent code
switch (machine)

    % for PCs/Windows hosts
    case {'WIN'}

        % for now, replace \'s with /'s
        tpath = strrep(tpath, '\', '/');
        % path can only be absolute for either of
        % - X:\<...>
        % - \\HOST\PATH_TO_UNC
        if (numel(tpath) > 2 && ...
            tpath(1) > 64 && ...
            tpath(1) < 123 && ...
            (tpath(1) < 91 || tpath(1) > 96) && ...
            strcmp(tpath(2:3), ':/')) || ...
           (numel(tpath) > 3 && ...
            strcmp(tpath(1:2), '//'))
            isabs = true;
        elseif tpath(1) == '/'
            pwdpath = pwd;
            abspath = [pwdpath(1) ':' tpath];
        else
            abspath = strrep([pwd filesep tpath], '//', '/');
        end

    % for Linux/Unix'ish systems
    case {'UNIX'}

        % path can only be "absolute" for either of
        % - /...
        % - ~/...
        if tpath(1) == '/' || ...
            (numel(tpath) > 1 && ...
             strcmp(tpath(1:2), '~/'))
            isabs = true;
        else
            abspath = strrep([pwd filesep tpath], '//', '/');
        end

    % for Macintosh
    case {'MAC'}

        % path could be absolute for either of
        % - /...
        % - ~/...
        if tpath(1) == '/' || ...
            (numel(tpath) > 1 && ...
             strcmp(tpath(1:2), '~/'))
            isabs = true;
        else
            abspath = strrep([pwd filesep tpath], '//', '/');
        end
end
