function [id, iw] = flexinterpnw_method(data, dw, cr, varargin)
% flexinterpnw_method  - call flexinterpnw with one of the typical methods
%
% FORMAT:       [id, iw] = flexinterpnw_method(d, dw, cr [, [opts], method])
%
% Input fields:
%
%       d           N-d data (up to 4 dims supported)
%       dw          N-d weights (double, same size as data)
%       cr          CxD coordinates or 4xD range
%       opts        hold (and tranformation matrix)
%       method      either of:
%                   'cubic', 'lanczos2' 'lanczos3', {'linear'},
%                   'nearest', 'poly3', 'spline2', 'spline3'
%
% Output fields:
%
%       id          interpolated data as from flexinterpn
%       iw          interpolated weights
%
% Note: the range notation of cr must contain Infs in the first row
%       of values, the rest is parsed as row2:row3:row4 for each dim.
%
% Note: for methods other than nearest and linear, the value range of
%       the input can be exceeded!

% Version:  v0.9c
% Build:    11120715
% Date:     Dec-07 2011, 3:36 PM EST
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

% Note: general information on kernel math taken from
% http://www.all-in-one.ee/~dersch/interpolator/interpolator.html

% compute kernels in persistent memory
persistent i_fwimkernel;
if isempty(i_fwimkernel)
    i_fwimkernel = flexinterpn_method('kernels');
end

% argument check
if nargin < 3 || ...
   ~isnumeric(data) || ...
   ~isa(dw, 'double') || ...
   ~isequal(size(data), size(dw)) || ...
   ~isnumeric(cr)
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid argument.' ...
    );
end
narg = nargin;
if narg < 4 || ...
   ~ischar(varargin{end}) || ...
    isempty(regexpi(varargin{end}(:)', '^(cubic|lanczos[2-9]|linear|nearest|poly3|spline[23])$'))
    method = 'linear';
else
    method = lower(varargin{end}(:)');
    narg = narg - 1;
end
if narg > 4 && ...
    ischar(varargin{end-1}) && ...
    strcmpi(varargin{end-1}(:)', 'gauss')
    if numel(varargin{end}) == 1 && ...
        isa(varargin{end}, 'double') && ...
       ~isnan(varargin{end}) && ...
        varargin{end} >= 0.25 && ...
        varargin{end} <= 8
        gaussk = varargin{end};
    else
        gaussk = 2;
    end
    method = 'gauss';
    narg = narg - 2;
end

% special case, nearest neighbor without rotation, up to 3D
if strcmp(method, 'nearest') && ...
    (narg < 4 || ...
     (narg == 4 && ...
      numel(varargin{1}) == 1)) && ...
    ndims(data) < 4

    % use indexarray
    try
        id = indexarraynb(data .* dw, cr);
        iw = indexarraynb(dw, cr);
        id = id ./ iw;
        id(iw == 0) = 0;
    catch ne_eo;
        rethrow(ne_eo);
    end
    return;
end

% depending on method
if ~strcmp(method, 'gauss')
    k = i_fwimkernel.(method);

% gaussian smoothing while interpolating
else
    f  = gaussk / sqrt(8 * log(2));
    md = round(6 * f);
    k = {exp(- (-md:1/1024:md)' .^ 2 ./ (2 * f .^ 2)), 1024};
    k{1}([1, numel(k{1})]) = 0;
end

% special case, rotation of 2D data
if ndims(data) == 2 && ...
    isequal(size(cr), [4, 2]) && ...
    all(isinf(cr(1, :))) && ...
    narg > 4 && ...
    isa(varargin{2}, 'double') && ...
    isequal(size(varargin{2}), [4, 4])

    % perform this in 3D!
    data = reshape(data, [1, size(data)]);
    cr = [[Inf;1;1;1], cr];

    % pass on and return
    try
        [id, iw] = shiftdim(flexinterpnw(data, dw, cr, k{:}, varargin{1:narg-3}), 1);
    catch ne_eo;
        rethrow(ne_eo);
    end
    return;
end

% now pass this on
try
    [id, iw] = flexinterpnw(data, dw, cr, k{:}, varargin{1:narg-3});
catch ne_eo;
    rethrow(ne_eo);
end
