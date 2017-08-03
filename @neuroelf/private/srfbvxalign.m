function succ = srfbvxalign(bvx)
% srfbvxalign  - perform alignment of SRFs in a BVX file
%
% FORMAT:       [succ = ] srfbvxalign(bvx)
%
% Input fields:
%
%       bvx         BVX file with necessary information
%
% Output fields:
%
%       succ        boolean success flag

% Version:  v1.1
% Build:    16012313
% Date:     Jan-23 2016, 1:55 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2016, Jochen Weber
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
if nargin ~= 1 || ...
   ~ischar(bvx) || ...
    numel(bvx) < 5 || ...
    numel(bvx) ~= size(bvx, 2) || ...
   ~strcmpi(bvx(end-3:end), '.bvx') || ...
    exist(bvx, 'file') ~= 2
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end

% load BVX with transio access
xffroot = xff();
btio = xffroot.TransIOSize('bvx', 4096);
obvx = xff(bvx);
xffroot.TransIOSize('bvx', btio);

% get options
nsrf = obvx.Variables(1).Content;
alst = obvx.Variables(2).Content;
icom = (obvx.Variables(3).Content ~= 0);
nvrt = obvx.Variables(4).Content;
ts_c = obvx.Variables(5).Content;
ts_v = obvx.Variables(6).Content(:, :);
ts_t = obvx.Variables(7).Content(:, :);
ts_m = obvx.Variables(8).Content(:, :);
if size(ts_t, 1) > 1
    [ts_n, ts_ne] = mesh_trianglestoneighbors(size(ts_v, 1), ts_t);
    if ~isempty(ts_ne)
        error( ...
            'neuroelf:InternalError', ...
            'Template surface seems to be corrupt.' ...
        );
    end
end

% all went well
succ = true;
