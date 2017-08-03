function cdata = scaledata(cdata, mm, mm2)
%SCALEDATA  Scale data within (grayscale) boundaries.
%   D = SCALEDATA(D) scales the ND array between 0 and 255.999.
%
%   D = SCALEDATA(D, BOUNDARIES) scales the ND array between BOUNDARIES(1)
%   and BOUNDARIES(2). For this to work, BOUNDARIES must be a 1x2 double
%   with two different values.
%
%   D = SCALEDATA(D, DIM) scales the data in D along DIM between 0 and 1.
%
%   D = SCALEDATA(D, DIM, BOUNDARIES) scales the data between boundaries.
%
%   See also PSCTRANS, ZTRANS.

% Version:  v1.0
% Build:    16032410
% Date:     Mar-24 2016, 10:29 AM EST
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

% argument check
if ~isnumeric(cdata)
    error('neuroelf:general:badArgument', 'Invalid color data given.');
end
cdata = double(cdata);

% dim given
dd = ndims(cdata);
if nargin > 1 && isa(mm, 'double') && numel(mm) == 1 && ~isnan(mm) && any(mm == (1:dd))

    % boundaries given as well
    if nargin < 3 || ~isa(mm2, 'double') || numel(mm2) ~= 2 || any(isinf(mm2) | isnan(mm2)) || mm2(1) == mm2(2)
        mm2 = [0, 1];
    end

    % compute scaling factor
    mmd = mm2(2) - mm2(1);

    % repmat argument
    rma = ones(1, dd);
    rma(mm) = size(cdata, mm);

    % remove min along dim
    cfac = min(cdata, [], mm);
    if numel(cdata) < 3.5e7
        cdata = cdata - repmat(cfac, rma);
    else
        sa = repmat({':'}, 1, dd);
        for dc = 1:size(cdata, mm)
            sa{mm} = dc;
            cdata(sa{:}) = cdata(sa{:}) - cfac;
        end
    end

    % multiply by correct scaling factor
    cfac = mmd ./ (eps + max(cdata, [], mm));
    if numel(cdata) < 3.5e7
        cdata = repmat(cfac, rma) .* cdata;
    else
        for dc = 1:size(cdata, mm)
            sa{mm} = dc;
            cdata(rma{:}) = cfac .* cdata(sa{:});
        end
    end

    % if necessary add offset
    if mm2(1) ~= 0
        cdata = mm2(1) + cdata;
    end

    % return
    return;
end

% regular process
if nargin < 2 || ~isa(mm, 'double') || numel(mm) ~= 2 || ...
    any(isinf(mm) | isnan(mm)) || mm(1) == mm(2)
    mm = [0, 255.999];
end

% compute factor
md = mm(2) - mm(1);

% and min/max
mmm = minmaxmean(cdata(:), 4);

% compute
cdata = mm(1) + (md / (eps + (mmm(2) - mmm(1)))) .* (cdata - mmm(1));
