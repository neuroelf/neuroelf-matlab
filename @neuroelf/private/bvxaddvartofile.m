function bvxaddvartofile(filename, varname, varcont)
% bvxaddvartofile  - add a variable to an existing BVX file
%
% FORMAT:       [succ = ] bvxaddvartofile(filename, varname, varcont)
%
% Input fields:
%
%       filename    BVX filename
%       varname     variable name(s), either 1xC char or 1xV cell array
%       varcont     variable content(s), either numeric or 1xV cell array
%
% No output fields.

% Version:  v0.9b
% Build:    11050712
% Date:     Apr-09 2011, 1:55 PM EST
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

% persistent type struct
persistent bvx_classid;
if isempty(bvx_classid)
    bvx_classid = struct( ...
        'double',  6, ...
        'single',  7, ...
        'int8',    8, ...
        'uint8',   9, ...
        'int16',  10, ...
        'uint16', 11, ...
        'int32',  12, ...
        'uint32', 13, ...
        'int64',  14, ...
        'uint64', 15);
end

% argument check
if nargin < 3 || ...
   ~ischar(filename) || ...
    numel(filename) < 5 || ...
    numel(filename) ~= size(filename, 2) || ...
   ~strcmpi(filename(end-3:end), '.bvx') || ...
   (~iscell(varname) && ...
    (~ischar(varname) || ...
     isempty(varname) || ...
     numel(varname) > 63 || ...
   ~strcmp(varname(:)', makelabel(varname(:)')))) || ...
   ((~iscell(varcont) || ...
     isempty(varcont)) && ...
   ~isnumeric(varcont)) || ...
   (iscell(varname) && ...
   ~iscell(varcont))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if ~iscell(varname)
    varname = {varname(:)'};
end
if ~iscell(varcont)
    varcont = {varcont};
end
if numel(varname) ~= numel(varcont)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end

% open file for append
fid = fopen(filename, 'a+', 'ieee-le');
if fid < 1
    error( ...
        'neuroelf:FileOpenError', ...
        'Error opening file for read/append IO.' ...
    );
end

% check minimum size
fseek(fid, 0, 1);
fl = ftell(fid);
fseek(fid, 0, -1);
if fl < 8
    fclose(fid);
    return;
end

% check fileversion
fv = fread(fid, [1, 1], 'uint32=>double');
if ~any(1 == fv)
    fclose(fid);
    return;
end

% read current number of variables
nv = fread(fid, [1, 1], 'uint32=>double');

% seek to end
if nv > 0
    fseek(fid, 0, 1);
end

% test variable classes
nnv = numel(varcont);
varcid = zeros(1, nnv);
for vc = 1:numel(varcont)
    try
        varcid(vc) = bvx_classid.(lower(class(varcont{vc})));
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        fclose(fid);
        error( ...
            'neuroelf:BadArgument', ...
            'Invalid class in varcont.' ...
        );
    end
end

% write variables
for vc = 1:numel(varcont)
    we = fwrite(fid, [double(varname{vc}(:)'), 0], 'uint8');
    if we ~= (numel(varname{vc}) + 1)
        fclose(fid);
        error( ...
            'neuroelf:FileIOError', ...
            'Error writing variable name of var %d.', ...
            vc ...
        );
    end
    we = fwrite(fid, [ndims(varcont{vc}), size(varcont{vc}), varcid([vc, vc])], 'uint32');
    if we ~= (ndims(varcont{vc}) + 3)
        fclose(fid);
        error( ...
            'neuroelf:FileIOError', ...
            'Error writing variable size/dims/class of var %d.', ...
            vc ...
        );
    end
    we = fwrite(fid, varcont{vc}, lower(class(varcont{vc})));
    if we ~= numel(varcont{vc})
        fclose(fid);
        error( ...
            'neuroelf:FileIOError', ...
            'Error writing variable content of var %d.', ...
            vc ...
        );
    end
end

% seek back to NrOfVariables
fclose(fid);
fid = fopen(filename, 'r+', 'ieee-le');
if fid < 1
    error( ...
        'neuroelf:FileOpenError', ...
        'Error re-opening file for updating number of variables.' ...
    );
end
fseek(fid, 4, -1);
if ftell(fid) ~= 4
    error( ...
        'neuroelf:FileOpenError', ...
        'Error seeking within file to rewrite number of variables.' ...
    );
end

% write new number of variables
we = fwrite(fid, nv + nnv, 'uint32');
fclose(fid);
if we ~= 1
    error( ...
        'neuroelf:FileIOError', ...
        'Error updating number of variables.' ...
    );
end
