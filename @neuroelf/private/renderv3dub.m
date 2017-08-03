function renderv3dub(data, vc)
% renderv3dub  - updating the internal render buffer (Ibuffer)
%
% FORMAT:       renderv3dub(data, vc)
%
% Input fields:
%
%       data        struct with required fields
%        .Ibuffer   image buffer (written into!)
%        .alpha_loc location-based alpha (if separate from intensity)
%        .atable    1xV cell array with alpha-table information
%        .ctable    1xV cell array with color-table information
%        .imaxmin   1xV scaling factor (max - min) per volume
%        .imin      1xV minimum value per volume
%        .intensity_loc  location-based intensities (unscaled)
%        .intensity_rgb  location-based RGB values (unscaled)
%        .px        X-positions in Ibuffer
%        .py        Y-positions in Ibuffer
%       vc          lookup index for volume (to access atable, ctable)
%
% No output fields.
%
% Note: for now this only supports RGB (color) updating!

% Version:  v0.9c
% Build:    14030510
% Date:     Mar-05 2014, 10:39 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2014, Jochen Weber
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

% error out without mex
error( ...
    'neuroelf:MissingMEXFile', ...
    'This is a compiled function, but the MEX file is missing.' ...
);
