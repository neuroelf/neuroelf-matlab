function succ = bvxcreatefile(filename)
% bvxcreatefile  - create new BVX file
%
% FORMAT:       [succ = ] bvxcreatefile(filename)
%
% Input fields:
%
%       filename    filename of BVX file to create
%
% Output fields:
%
%       succ        success

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:14 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2016, Jochen Weber
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

% default
succ = false;

% argument check
if nargin < 1 || ...
   ~ischar(filename) || ...
    numel(filename) < 5 || ...
    numel(filename) ~= size(filename, 2)
    return;
end
if ~strcmpi(filename(end-3:end), '.bvx')
    return;
end

% file already exists
if exist(filename, 'file') == 2
    try
        delete(filename);
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
    end
    if exist(filename, 'file') == 2
        return;
    end
end

% try creation
bvx = [];
try
    bvx = xff('new:bvx');
    bvx.SaveAs(filename);
catch ne_eo;
    if isxff(bvx, true)
        bvx.ClearObject;
    end
    rethrow(ne_eo);
end
if ~isempty(bvx)
    try
        bvx.ClearObject;
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
    end
end
if exist(filename, 'file') == 2
    succ = true;
end
