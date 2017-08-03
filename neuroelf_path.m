function p = neuroelf_path(p)
% neuroelf_path  - get path where NeuroElf is installed
%
% FORMAT:       p = neuroelf_path([p])
%
% Input fields:
%
%       p           optional path, either of
%                   - base (default)
%                   - bin
%                   - cache
%                   - colin
%                   - config
%                   - core
%                   - fonts
%                   - formats
%                   - files
%                   - icons
%                   - images
%                   - lut
%                   - masks
%                   - pines
%                   - remote
%                   - shen
%                   - splash
%                   - spm
%                   - srf
%                   - tal
%                   - tfg
%
% Output fields:
%
%       p           path with forward slashes for xff (or other)

% Version:  v1.1
% Build:    16060810
% Date:     Jun-08 2016, 10:27 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2016, Jochen Weber
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

% persistent memory
persistent nelf_p;
if numel(nelf_p) ~= 1 || ~isstruct(nelf_p)

    % create struct
    nelf_p = struct;

    % get path
    px = strrep(fileparts(mfilename('fullpath')), '\', '/');
    if px(end) == '/'
        px(end) = [];
    end

    % put path in base field
    nelf_p.base = px;
    p2 = [px '/_core'];
    p3 = [px '/_files'];

    % other folders
    nelf_p.bin     = [p3 '/compiled/' strrep(mexext, 'mex', '')];
    nelf_p.cache   = [p3 '/cache'];
    nelf_p.colin   = [p3 '/colin'];
    nelf_p.config  = [p2 '/config'];
    nelf_p.core    = p2;
    nelf_p.dicom   = [p2 '/dicom'];
    nelf_p.fonts   = [p2 '/fonts'];
    nelf_p.formats = [p2 '/formats'];
    nelf_p.files   = p3;
    nelf_p.icons   = [p2 '/icons'];
    nelf_p.images  = [p3 '/images'];
    nelf_p.lut     = [p2 '/lut'];
    nelf_p.masks   = [p3 '/masks'];
    nelf_p.nsynth  = [p3 '/neurosynth'];
    nelf_p.pines   = [p3 '/pines'];
    nelf_p.remote  = [p3 '/contrib/ne_remote'];
    nelf_p.shen    = [p3 '/shenparcel'];
    nelf_p.splash  = [p2 '/splash'];
    nelf_p.spm     = [p3 '/spm'];
    nelf_p.srf     = [p3 '/srf'];
    nelf_p.tal     = [p3 '/tal'];
    nelf_p.tfg     = [p2 '/tfg'];
end

% argument check
if nargin < 1
    p = 'base';
elseif ~ischar(p) || isempty(p)
    error('neuroelf:general:badArgument', 'Invalid p argument.');
else
    p = lower(p(:)');
end

% pre-set
if isfield(nelf_p, p)
    p = nelf_p.(p);
    return;
end

% found in files
if exist([nelf_p.files '/' p], 'dir') > 0
    p = [nelf_p.files '/' p];

% else return base path
else
    warning('neuroelf:general:badArgument', ...
        'Invalid pathspec: %s, returning base dir.', p);
    p = nelf_p.base;
end
