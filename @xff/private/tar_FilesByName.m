function [f, p, l, t] = tar_FilesByName(xo, n, opts)
% TAR::FilesByName  - look for files by name
%
% FORMAT:       [files, pos, lengths, tio] = tar.FilesByName(name [, opts])
%
% Input fields:
%
%       name        string (1xN char) identifying which files to select
%       opts        1x1 struct with optional settings
%        .icase     1x1 logical flag, ignore case (default: true)
%        .transio   1x1 logical flag, return cell of transios (false)
%
% Output fields:
%
%       files       filenames that match
%       pos         position within the TAR file
%       lengths     file lengths
%       tio         cell array with transio objects

% Version:  v1.1
% Build:    16060314
% Date:     Jun-03 2016, 2:56 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
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

% only valid for single file
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'tar') || ...
   ~ischar(n) || isempty(n)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
n = n(:)';
if ~any(n == '*') && ~any(n == '?') && ~any(n == '|')
    n = [n '$'];
end

% options
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'icase') || ~islogical(opts.icase) || numel(opts.icase) ~= 1
    opts.icase = true;
end
if ~isfield(opts, 'transio') || ~islogical(opts.transio) || numel(opts.transio) ~= 1
    opts.transio = false;
end

% match
if opts.icase
    m = find(~cellfun('isempty', regexpi(xo.C.Name, n)));
else
    m = find(~cellfun('isempty', regexp(xo.C.Name, n)));
end

% nothing found
if isempty(m)
    f = cell(0, 1);
    p = zeros(0, 1);
    l = zeros(0, 1);
    t = cell(0, 1);

% found
else

    % remove folders
    m(xo.C.Type(m) > '0') = [];

    % get files, positions and lengths
    f = xo.C.Name(m);
    p = xo.C.ContPos(m);
    l = xo.C.ContLen(m);
    if nargout > 3
        t = cell(numel(m), 1);

        % transios?
        if opts.transio
            tio_obj = struct(transio(xo.F, 'ieee-le', 'uint8', 0, [1, 512]));
            tio_pos = p - 1;
            for fc = 1:numel(m)
                tio_obj.IOOffset = tio_pos(fc);
                tio_obj.DataDims = [1, l(fc)];
                t{fc} = transio(0, 'makeobject', tio_obj);
            end
        end
    end
end

