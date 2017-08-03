function [h, s] = hrf(shape, sfrq, pttp, nttp, pnr, ons, pdsp, ndsp, s, d, nt)
% hrf  - create canonical HRF shape
%
% FORMAT:       [h, s] = hrf(shape, sf, pttp, nttp, pnr, ons, pdsp, ndsp, s, d, nt)
%
% Input fields:
%
%       shape       HRF general shape {'twogamma' [, 'boynton']}
%       sf          HRF sample frequency (default: 1s/16, OK: [1e-3 .. 5])
%       pttp        time to positive (response) peak (default: 5 secs)
%       nttp        time to negative (undershoot) peak (default: 15 secs)
%       pnr         pos-to-neg ratio (default: 6, OK: [1 .. Inf])
%       ons         onset of the HRF (default: 0 secs, OK: [-5 .. 5])
%       pdsp        dispersion of positive gamma PDF (default: 1)
%       ndsp        dispersion of negative gamma PDF (default: 1)
%       s           sampling range (default: [0, ons + 2 * (nttp + 1)])
%       d           derivatives (default: 0)
%       nt          either of {'area'} or 'sum'
%
% Output fields:
%
%       h           HRF function given within [0 .. onset + 2*nttp]
%       s           HRF sample points
%
% Note: the pttp and nttp parameters are increased by 1 before given
%       as parameters into the gammapdf function (which is a property
%       of the gamma PDF!)
%
% Note: SPM normalizes the functions based on their sum, whereas
%       the default for NeuroElf is to normalize by area (so as to allows
%       the Calhoun-based HRF-boosting!)
%
% See also gammapdf

% Version:  v0.9c
% Build:    12041921
% Date:     Apr-12 2012, 5:19 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2012, Jochen Weber
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

% argument check
if nargin < 10 || ...
   ~isa(d, 'double') || ...
    numel(d) ~= 1 || ...
    isinf(d) || ...
    isnan(d) || ...
    d < 0 || ...
    d > 2
    d = 0;
else
    d = round(d);
end
if nargin < 8 || ...
   ~isa(ndsp, 'double') || ...
    numel(ndsp) ~= 1 || ...
    isinf(ndsp) || ...
    isnan(ndsp) || ...
    ndsp <= 0
    ndsp = 1;
end
if nargin < 7 || ...
   ~isa(pdsp, 'double') || ...
    numel(pdsp) ~= 1 || ...
    isinf(pdsp) || ...
    isnan(pdsp) || ...
    pdsp <= 0
    pdsp = 1;
end
if nargin < 6 || ...
   ~isa(ons, 'double') || ...
    numel(ons) ~= 1 || ...
    isinf(ons) || ...
    isnan(ons) || ...
    ons < -5 || ...
    ons > 5
    ons = 0;
end
if nargin < 5 || ...
   ~isa(pnr, 'double') || ...
    numel(pnr) ~= 1 || ...
    isnan(pnr) || ...
    pnr < 1
    pnr = 6;
end
if nargin < 4 || ...
   ~isa(nttp, 'double') || ...
    numel(nttp) ~= 1 || ...
    isinf(nttp) || ...
    isnan(nttp) || ...
    nttp < 1 || ...
    nttp > 30
    nttp = 15;
end
if nargin < 3 || ...
   ~isa(pttp, 'double') || ...
    numel(pttp) ~= 1 || ...
    isinf(pttp) || ...
    isnan(pttp) || ...
    pttp < 1 || ...
    pttp > 30
    pttp = 5;
end
if nargin < 2 || ...
   ~isa(sfrq, 'double') || ...
    numel(sfrq) ~= 1 || ...
    isinf(sfrq) || ...
    isnan(sfrq) || ...
    sfrq > 5
    sfrq = 1/16;
elseif sfrq < 0.001
    sfrq = 0.001;
end
if nargin < 1 || ...
   ~ischar(shape) || ...
   ~any(strcmpi(shape, {'twogamma', 'boynton'}))
    shape = 'twogamma';
else
    shape = lower(shape(:)');
end
if nargin < 9 || ...
   ~isa(s, 'double') || ...
    numel(s) ~= 2 || ...
    any(isinf(s) | isnan(s))
    s = (0:sfrq:(ons + 2 * (nttp + 1))) - ons;
else
    s = min(s):sfrq:max(s) - ons;
end
if nargin < 11 || ...
   ~ischar(nt) || ...
    isempty(nt) || ...
    lower(nt(1)) ~= 's'
    nt = 'a';
else
    nt = lower(nt(1));
end

% computation (according to shape)
h = zeros(numel(s), d + 1);
switch (shape)

    % boynton (single-gamma) HRF
    case {'boynton'}
        h(:, 1) = gammapdf(s(:), pttp + 1, pdsp);
        if d > 0
            h(:, 2) = h(:, 1) - gammapdf(s(:) + 1, pttp + 1, pdsp);
            hi = find(h(:, 2) ~= 0);
            h(hi, 2) = h(hi, 2) - ((pinv(h(hi, 1)' * h(hi, 1)) * h(hi, 1)' * h(hi, 2))' * h(hi, 1)')';
            if d > 1
                h(:, 3) = h(:, 1) - gammapdf(s(:), pttp + 1, pdsp / 1.01);
                hi = find(h(:, 3) ~= 0);
                h(hi, 3) = h(hi, 3) - ((pinv(h(hi, [1, 2])' * h(hi, [1, 2])) * ...
                    h(hi, [1, 2])' * h(hi, 3))' * h(hi, [1, 2])')';
            end
        end

    % two-gamma HRF
    case {'twogamma'}
        h(:, 1) = gammapdf(s(:), pttp + 1, pdsp) - ...
            gammapdf(s(:), nttp + 1, ndsp) / pnr;
        if d > 0
            h(:, 2) = h(:, 1) - (gammapdf(s(:) - 1, pttp + 1, pdsp) - ...
                gammapdf(s(:) - 1, nttp + 1, ndsp) / pnr);
            hi = find(h(:, 2) ~= 0);
            h(hi, 2) = h(hi, 2) - ((pinv(h(hi, 1)' * h(hi, 1)) * h(hi, 1)' * h(hi, 2))' * h(hi, 1)')';
            if d > 1
                h(:, 3) = (h(:, 1) - (gammapdf(s(:), (pttp + 1) / 1.01, pdsp / 1.01) - ...
                    gammapdf(s(:), nttp + 1, ndsp) / pnr)) ./ 0.01;
                hi = find(h(:, 3) ~= 0);
                h(hi, 3) = h(hi, 3) - ((pinv(h(hi, [1, 2])' * h(hi, [1, 2])) * ...
                    h(hi, [1, 2])' * h(hi, 3))' * h(hi, [1, 2])')';
            end
        end
end

% normalize for area first
if nt == 'a' && ...
    d > 0
    h = h ./ (ones(size(h, 1), 1) * sqrt(sum(h .* h)));
end

% then normalize everything
h = h ./ sum(h(:, 1));
