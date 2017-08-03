function talout = icbm2tal(icbmin, rounded)
% icbm2tal - converts coordinates from segmented ICBM space to TAL space
%
% FORMAT:       talout = icbm2tal(icbmin [, rounded])
%
% Input fields:
%
%       icbmin      N-by-3 or 3-by-N matrix of coordinates
%       rounded     1x1 double, if given, coordinates are rounded to
%                   specified number of digits
%
% Output fields:
%
%       talout      is the coordinate matrix with TAL-space points

% Version:  v0.9b
% Build:    11050712
% Date:     Apr-09 2011, 11:08 PM EST
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

% persistent transformation matrices
persistent i2t_trf;
if isempty(i2t_trf) || ...
   ~isstruct(i2t_trf)
    try
        i2t_trf = load(neuroelf_file('t', 'talairach_seg_sn.mat'));
        if ~isfield(i2t_trf, 'Affine') || ...
           ~isfield(i2t_trf, 'VG') || ...
           ~isstruct(i2t_trf.VG) || ...
            isempty(i2t_trf.VG) || ...
           ~isfield(i2t_trf.VG, 'dim') || ...
           ~isfield(i2t_trf.VG, 'mat')
            error( ...
                'neuroelf:BadFileContent', ...
                'Talairach->ICBM normalization file invalid or missing.' ...
            );
        end
    catch ne_eo;
        i2t_trf = [];
        rethrow(ne_eo);
    end
    i2t_trf.VG = i2t_trf.VG(1);
    i2t_trf.VGmat = inv(i2t_trf.VG.mat);
end

% argument check
if nargin < 1 || ...
    length(size(icbmin)) > 2 || ...
   ~isa(icbmin, 'double')
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing input argument mniin.' ...
    );
end

% transpose ?
dimdim = find(size(icbmin) == 3);
if isempty(dimdim)
    error( ...
        'neuroelf:BadArguments', ...
        'talin argument must be a N-by-3 or 3-by-N matrix' ...
    );
end

% transpose as needed
if dimdim(1) == 1 && ...
    numel(dimdim) == 1
    icbmin = icbmin';
end

% use applyspmsnc for transformation
talout = applyspmsnc(icbmin, i2t_trf.Tr, i2t_trf.VG.dim, ...
    i2t_trf.VGmat, i2t_trf.VF(1).mat * i2t_trf.Affine);

% retranspose to match input
if dimdim(1) == 1 && ...
    numel(dimdim) == 1
    talout = talout';
end

% round output
if nargin > 1 && ...
    isa(rounded, 'double') && ...
   ~isempty(rounded)
    talout = (1/10^rounded(1)) * round(talout * 10^rounded(1));
end
