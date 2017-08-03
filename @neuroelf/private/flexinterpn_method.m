function [idata, k] = flexinterpn_method(data, cr, varargin)
% flexinterpn_method  - call flexinterpn with one of the typical methods
%
% FORMAT:       [idata, k] = flexinterpn_method(data, cr [, [opts], method])
%
% Input fields:
%
%       data        N-d data (up to 4 dims supported)
%       cr          CxD coordinates or 4xD range
%       opts        hold (and tranformation matrix)
%       method      either of:
%                   'cubic', 'lanczos2' 'lanczos3', {'linear'},
%                   'nearest', 'poly3', 'spline2', 'spline3'
%
% Output fields:
%
%       idata       interpolated data as from flexinterpn
%       k           kernel passed to flexinterpn (empty for linear)
%
% Note: the range notation of cr must contain Infs in the first row
%       of values, the rest is parsed as row2:row3:row4 for each dim.
%
% Note: for methods other than nearest and linear, the value range of
%       the input can be exceeded!

% Version:  v0.9b
% Build:    14011714
% Date:     Apr-09 2011, 2:11 PM EST
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
persistent i_fimkernel;
if isempty(i_fimkernel)
    i_fimkernel = struct;
    k21 = 4097:12289;
    k22 = [1:4096, 12290:16385];

    % cubic spline, see http://en.wikipedia.org/wiki/Bicubic_interpolation
    % create line from -2 ... 2 (with 1/4096 steps)
    k = abs(-2:1/4096:2)';

    % compute term for 0 <= |x| <= 1
    k(k21) =  1.5 * (k(k21) .^ 3) - 2.5 * (k(k21) .^ 2) + 1.0;

    % compute term for 1 < |x| <= 2
    k(k22) = -0.5 * (k(k22) .^ 3) + 2.5 * (k(k22) .^ 2) - 4.0 * k(k22) + 2.0;

    % set all "direct hit" nodes
    k(1:4096:end) = [0, 0, 1, 0, 0];

    % store in persistent struct
    i_fimkernel.cubic = {k, 4096};


    % Lanczos kernels, see http://en.wikipedia.org/wiki/Lanczos_resampling
    % from kernel size 2:9
    kss = [0, 8192, 8192, 4096, 4096, 2048, 2048, 2048, 2048];
    for kc = 2:9

        % start with line again
        ksi = kss(kc);
        k = (-kc:1/ksi:kc)';

        % set 0-value to 1 first
        k(kc * ksi + 1) = 1;

        % compute sinc / (pi * k)
        ks = sin(pi * k) ./ ((pi * k) .^ 2);
        ks(1:ksi:end) = 0;

        % Lanczos' addition
        ka = (kc * ks) .* sin((pi / kc) * k);

        % reset 0-value
        ka(kc * ksi + 1) = 1;

        % and store it
        i_fimkernel.(sprintf('lanczos%d', kc)) = {ka, ksi};
    end

    % store linear kernel as two empty arguments
    i_fimkernel.linear = {[], []};

    % for nearest neighbor, use pseudo sampling
    i_fimkernel.nearest = {[zeros(2048, 1); ones(4097, 1); zeros(2048, 1)], 4096};

    % poly3
    k = abs(-2:1/4096:2)';
    k(k21) = 1 + (-2.25 + 1.25 * k(k21)) .* (k(k21) .* k(k21));
    k(k22) = 3 + ((3.75 - 0.75 * k(k22)) .* (k(k22)) - 6) .* k(k22);
    k(1:4096:end) = [0,0,1,0,0]';
    i_fimkernel.poly3 = {k, 4096};

    % spline2 (diff not contiguous, but acceptable approx.)
    k = abs(-2:1/4096:2)';
    k(k21) = 1 + ((-1.8 + k(k21)) .* k(k21) - 0.2) .* k(k21);
    k(k22) = (((-1/3) * (-1 + k(k22)) + 0.8) .* (-1 + k(k22)) - (7/15)) .* (-1 + k(k22));
    k(1:4096:end) = [0,0,1,0,0]';
    i_fimkernel.spline2 = {k, 4096};

    % spline3 (diff not contiguous, but good approx.)
    k = abs(-3:1/4096:3)';
    k31 = 8193:16385;
    k32 = [4097:8192, 16386:20481];
    k33 = [1:4096, 20482:24577];
    k(k31) = 1 + (((13/11) * k(k31) - (453/209)) .* k(k31) - (3/209)) .* k(k31);
    k(k32) = (((-6/11) * (-1 + k(k32)) + (270/209)) .* (-1 + k(k32)) - (156/209)) .* (-1 + k(k32));
    k(k33) = (((1/11) * (-2 + k(k33)) - (45/209)) .* (-2 + k(k33)) + (26/209)) .* (-2 + k(k33));
    k(1:4096:end) = [0, 0, 0, 1, 0, 0, 0]';
    i_fimkernel.spline3 = {k, 4096};

end

% argument check
if nargin < 2 || ...
   ~isnumeric(data) || ...
   ~isnumeric(cr)

    % allow to access all kernels
    if nargin == 1 && ...
        ischar(data) && ...
       ~isempty(data) && ...
        strcmpi(data(:)', 'kernels')
        idata = i_fimkernel;
        return;
    end
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid argument.' ...
    );
end
narg = nargin;
if narg < 3 || ...
   ~ischar(varargin{end}) || ...
    isempty(regexpi(varargin{end}(:)', '^(cubic|lanczos[2-9]|linear|nearest|poly3|spline[23])$'))
    method = 'linear';
else
    method = lower(varargin{end}(:)');
    narg = narg - 1;
end
if narg > 3 && ...
    ischar(varargin{end-1}) && ...
    strcmpi(varargin{end-1}(:)', 'gauss')
    if numel(varargin{end}) == 1 && ...
        isa(varargin{end}, 'double') && ...
       ~isnan(varargin{end}) && ...
        varargin{end} >= 0.25 && ...
        varargin{end} <= 16
        gaussk = varargin{end};
    else
        gaussk = 2;
    end
    method = 'gauss';
    narg = narg - 2;
end

% special case, nearest neighbor without rotation, up to 3D
if strcmp(method, 'nearest') && ...
    (narg < 3 || ...
     (narg == 3 && ...
      numel(varargin{1}) == 1)) && ...
    ndims(data) < 4

    % use indexarray
    try
        if isequal(size(cr), [4, 3]) && ...
            all(isinf(cr(1, :)))
            cr(4, :) = cr(4, :) + 0.1 .* cr(3, :);
        end
        idata = indexarraynb(data, cr);
    catch ne_eo;
        rethrow(ne_eo);
    end
    return;
end

% depending on method
if ~strcmp(method, 'gauss')
    k = i_fimkernel.(method);

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
    narg > 3 && ...
    isa(varargin{2}, 'double') && ...
    isequal(size(varargin{2}), [4, 4])

    % perform this in 3D!
    data = reshape(data, [1, size(data)]);
    cr = [[Inf;1;1;1], cr];

    % pass on and return
    try
        idata = shiftdim(flexinterpn(data, cr, k{:}, varargin{1:narg-2}), 1);
    catch ne_eo;
        rethrow(ne_eo);
    end
    return;
end

% now pass this on
try
    idata = flexinterpn(data, cr, k{:}, varargin{1:narg-2});
catch ne_eo;
    rethrow(ne_eo);
end
