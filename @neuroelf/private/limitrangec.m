function V = limitrangec(V, rmin, rmax, rinv)
% limitrangec  - limit data to a range
%
% FORMAT:       V = limitrangec(V, rmin, rmax [, rinv])
%
% Input fields:
%
%       V           N-D numeric matrix
%       rmin        minimum value
%       rmax        maximum value
%       rinv        if given, replace invalid values (Inf/NaN) with this
%
% Output fields:
%
%       V           limited output
%
% Note: this is a MEX (compiled) function!

% Version:  v0.9c
% Build:    11042919
% Date:     Apr-29 2011, 8:11 PM EST
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

% alternative implementation
if nargin < 3 || ...
   ~isnumeric(V) || ...
   ~isnumeric(rmin) || ...
    numel(rmin) ~= 1 || ...
   ~isnumeric(rmax) || ...
    numel(rmax) ~= 1
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end

% only perform extended check for single/double
if nargin > 3 && ...
   (isa(V, 'double') || ...
    isa(V, 'single')) && ...
    isnumeric(rinv) && ...
    numel(rinv) == 1
    V(isinf(V) | isnan(V)) = rinv;
end

% simply use min/max
V(V < rmin) = rmin;
V(V > rmax) = rmax;
