function c = mesh_morph(c, n, tri, opts)
% mesh_morph -  morph the coordinates of a mesh
%
% FORMAT:       c = mesh_morph(c, n, tri, opts)
%
% Input fields:
%
%       c           Cx3 coordinate list (double)
%       n           Cx2 neighbors list (from SRF object, 1-based !)
%       tri         Tx3 triangle list (1-based !)
%       opts        mandatory struct with settings
%        .force     morphing force (1x1 double)
%        .niter     number of iterations (1x1 double)
%                 - optionally provided settings
%        .areac     if 1x1 double := 1, keep area constant
%                   (from initial state, requires .tri to be set!)
%        .distc     if given and between 0 .. 1, perform distortion corr
%        .distw     if given and [1], perform smoothing with distance
%                   weighting (default: false)
%        .distwsq   weight by square of distance (false)
%        .sphere    to-sphere force
%        .type      1xN char type, currently only 'smooth' supported
%
% Output fields:
%
%       c           morphed coordinates
%
% This is a MEX (c compiled) function for efficiency.

% Version:  v0.9d
% Build:    14062412
% Date:     Jun-24 2014, 12:50 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2014, Jochen Weber
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

% just to make sure
error( ...
    'neuroelf:MEXMissing', ...
    'This is a compiled function, but the MEX file is missing.' ...
);
