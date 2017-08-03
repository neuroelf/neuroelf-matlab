function [varargout] = icc(M, fz, ci, ft, fp)
% icc  - intraclass correlation calculation
%
% FORMAT:       [iccorr, ...] = icc(M [, fz, ci, ft])
%
% Input fields:
%
%       M           NxSxK array, where N = number of observations,
%                   S, number of subjects in group (ie. 2 for a pair),
%                   and K for multiple correlation calculations
%       fz          if given and true, apply Fisher's z-transform
%       ci          if given and true, also give confidence intervals
%       ft          if given and true, also give F-statistic
%       fp          if given and true, also give probabilities
%
% Output fields:
%
%       iccorr      Kx1 vector with intraclass correlations

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
   ~isa(M, 'double') || ...
    isempty(M) || ...
   ~isreal(M)
    error( ...
        'neuroelf:BadArgument', ...
        'Missing or bad input argument.' ...
    );
end

% sizes
size_M1 = size(M, 1);
size_M2 = size(M, 2);

% means
groupmean_M = mean(M, 2);
totalmean_M = mean(groupmean_M, 1);

% ANOVA calculation
total_ss__M = sum(sum( (M - totalmean_M(ones(1, size_M1), ones(1, size_M2), :)) .^ 2, 2), 1);
withinss__M = sum(sum( (M - groupmean_M(     1: size_M1 , ones(1, size_M2), :)) .^ 2, 2), 1);

betws_ss__M = total_ss__M - withinss__M;

mean_wss__M = withinss__M / (size_M1 * (size_M2 - 1));
mean_bss__M = betws_ss__M / (size_M1 - 1);

% squeeze result
iccorr = squeeze((mean_bss__M - mean_wss__M) ./ (mean_bss__M + mean_wss__M));

% transformations ?
if nargin < 2 || ...
    isempty(fz) || ...
   ~fz(1)
    fz = false;
else
    fz = true;
end
if nargin < 3 || ...
    isempty(ci) || ...
   ~ci(1)
    ci = false;
else
    ci = true;
end
if nargin < 4 || ...
    isempty(ft) || ...
   ~ft(1)
    ft = false;
else
    ft = true;
end
if nargin < 5 || ...
    isempty(fp) || ...
   ~fp(1)
    fp = false;
else
    fp = true;
end

% what exactly
if ft || fp
    ift = squeeze(mean_bss__M ./ mean_wss__M);
    if fp
        ifp = 1 - sdist('fcdf', ift, size_M1 - 1, size_M1 * (size_M2 - 1));
    end
end
if fz || ...
   (ci && size_M1 > 3)
    icc_z = 0.5 * log( (1 + iccorr) ./ (1 - iccorr) );
    if ci
        icc_zlow = icc_z - (1.96 / sqrt(size_M1 - 3));
        icc_zupp = icc_z + (1.96 / sqrt(size_M1 - 3));
    end
    if ~fz && ...
        ci
        e = exp(1);
        icc_zlow = ( (e .^ (2 * icc_zlow) - 1) ./ (e .^ (2 * icc_zlow) + 1) );
        icc_zupp = ( (e .^ (2 * icc_zupp) - 1) ./ (e .^ (2 * icc_zupp) + 1) );
    elseif fz
        iccorr = icc_z;
    end
end

% more output?
varargout{1} = iccorr;
if ci
    varargout{end+1} = icc_zlow;
    varargout{end+1} = icc_zupp;
end
if ft
    varargout{end+1} = ift;
end
if fp
    varargout{end+1} = ifp;
end
