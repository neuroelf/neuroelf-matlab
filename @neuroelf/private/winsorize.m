function [data, outl] = winsorize(data, sd, p)
% winsorize  - Winsorize data (to numbers of SD by numbers of passes)
%
% FORMAT:       [wd, outl] = winsorize(data [, sd [, p]])
%
% Input fields:
%
%       data        N-D data, along last dimension
%       sd          number of standard deviations (default: 3)
%       p           number of passes (default: 3)
%
% Output fields:
%
%       wd          winsorized data
%       outl        outliers (boolean array)

% Version:  v0.9b
% Build:    10081811
% Date:     Jun-17 2010, 12:06 AM EST
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
   ~isnumeric(data) || ...
    isempty(data) || ...
    any(isinf(data(:)) | isnan(data(:))) || ...
   ~isreal(data)
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid argument supplied.' ...
    );
end
if numel(data) == size(data, 1)
    mdim = 1;
else
    mdim = ndims(data);
end
if size(data, mdim) < 4
    return;
end
szd = size(data);
sdm = szd(mdim);
szd(mdim) = 1;
sdi = 1 / sdm;
if nargin < 2 || ...
    isempty(sd)
    sd = 3;
elseif ~isa(sd, 'double') || ...
    numel(sd) ~= 1 || ...
    isinf(sd) || ...
    isnan(sd) || ...
    sd <= 0.5
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid sd argument supplied.' ...
    );
end
if nargin < 3 || ...
    isempty(p)
    p = 3;
elseif ~isa(p, 'double') || ...
    numel(p) ~= 1 || ...
    isinf(p) || ...
    isnan(p) || ...
    p < 0 || ...
    p ~= fix(p)
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid p argument supplied.' ...
    );
end

% make sure data is double
dclass = class(data(1));
if ~isa(data, 'double')
    data = double(data);
end

% create boolean outlier argument
outl = false(size(data));
outs = (sdm - 1) .* ones(szd);

% create repmat argument
rma = ones(1, max(2, mdim));
rma(mdim) = sdm;

% compute mean
mdata = repmat(sdi .* sum(data, mdim), rma);

% iterate p times
for pc = 1:p

    % compute biased var
    sdata = varc(data, mdim, 1);

    % and then unbiased! std
    sdata = sqrt(((sdm - 1) ./ outs) .* sdata);

    % find entries that are further away than sd numbers of sd
    outl = (outl | (abs(data - mdata) > (sd .* repmat(sdata, rma))));

    % set outliers to mean
    data(outl) = mdata(outl);

    % recompute outs
    outsn = (sdm - 1) - sum(outl, mdim);
    if all(outsn(:) == outs(:))
        break;
    else
        outs = outsn;
    end

    % recompute mean
    mdata = repmat(sdi .* sum(data, mdim), rma);
end

% put back into old class?
if ~strcmpi(dclass, 'double')
    switch(dclass)
        case {'single'}
            data = single(data);
        case {'int32'}
            data = int32(round(data));
        case {'int16'}
            data = int16(round(data));
        case {'int8'}
            data = int8(round(data));
        case {'int64'}
            data = int64(round(data));
        case {'uint32'}
            data = uint32(round(data));
        case {'uint16'}
            data = uint16(round(data));
        case {'uint8'}
            data = uint8(round(data));
        case {'uint64'}
            data = uint64(round(data));
        case {'logical'}
            data = logical(data);
    end
end
