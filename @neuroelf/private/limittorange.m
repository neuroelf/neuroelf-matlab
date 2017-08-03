function [ld, minv, maxv] = limittorange(data, range, torange)
% limittorange  - limit data to percentile range
%
% FORMAT:       [ld, minv, maxv] = limittorange(data, range)
%
% Input fields:
%
%       data        numeric data
%       range       percentile (e.g. [0.05, 0.95]) or SDs (default: 3)
%       torange     if given also 1x2 numerical map on to numeric range
%
% Output fields:
%
%       ld          limited data
%       minv        minimum value
%       maxv        maximum value

% Version:  v0.9a
% Build:    10051716
% Date:     May-17 2010, 10:48 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, Jochen Weber
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
if nargin < 1 || ...
   ~isnumeric(data)
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid or missing argument in call.' ...
    );
end
if isempty(data)
    ld = data;
    return;
end
if nargin < 2 || ...
   ~isnumeric(range) || ...
   ~any(numel(range) == [1, 2]) || ...
    any(isinf(range) | isnan(range))
    range = 3;
else
    range = abs(double(range));
end
if numel(range) == 1 && ...
    numel(data) < 2
    ld = data;
    return;
else
    range = max(min(sort(range(:)), 1), 0);
    if all(range == [0; 1])
        ld = data;
        return;
    end
end

ld = data;
if numel(range) == 2
    if ~isa(data, 'double') && ...
       ~isa(data, 'single')
        data = single(data);
    end
    mmm = minmaxmean(ld, 4);
    mmh = (mmm(2) - mmm(1)) / 499;
    range = range .* numel(data);
    [h, hb] = hist(data(:), mmm(1):mmh:mmm(2));
    hs = cumsum(h);
    minh = findfirst(hs >= range(1));
    maxh = findfirst(hs <= range(2), -1);
    minv = data(1);
    maxv = minv;
    minv(1) = hb(minh);
    maxv(1) = hb(maxh);
else
    mmm = minmaxmean(data, 5);
    msd = sqrt(mmm(6));
    minv = data(1);
    maxv = minv;
    minv(1) = mmm(3) - range * msd;
    maxv(1) = mmm(3) + range * msd;
end
ld = min(max(ld, minv), maxv);
if nargin > 2 && ...
    isnumeric(torange) && ...
    numel(torange) == 2 && ...
   ~any(isinf(rorange) | isnan(torange))
    torange = double(torange);
    mmm = minmaxmean(ld, 4);
    ld = torange(1) + ((torange(2) - torange(1)) ./ (mmm(2) - mmm(1))) .* ...
        double(ld - mmm(1));
end
