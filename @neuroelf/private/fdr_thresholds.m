function th = fdr_thresholds(plist, levels, lnmeth, fdrfac)
% fdr_thresholds  - compute probability thresholds for q(FDR) levels
%
% FORMAT:       th = fdr_thresholds(plist, levels [, lnmeth [, fdrfac]])
%
% Input fields:
%
%       plist       list of p-values to be thresholded
%       levels      q(FDR) levels
%       lnmeth      method to use, default: 3
%                   0 :  c(V) = 1
%                   1 :  c(V) = ln(V) + E
%                   2 :  c(V) = 1, use last crossing
%                   3 :  c(V) = ln(V) + E, use last crossing
%       fdrfac      FDR factor (default: number of values)
%
% Output fields:
%
%       th          p-thresholds for q(FDR)

% Version:  v0.9c
% Build:    12011212
% Date:     Jan-05 2012, 11:29 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2012, Jochen Weber
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
   (~isa(plist, 'double') && ...
    ~isa(plist, 'single')) || ...
    numel(plist) ~= length(plist) || ...
    any(plist < 0)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing first argument.' ...
    );
end

% levels
if nargin < 2 || ...
   (~isa(levels, 'double') && ...
    ~isa(levels, 'single')) || ...
    isempty(levels) || ...
    numel(levels) ~= length(levels) || ...
    any(levels < 0 | levels >= 1)

    % assume 0.05
    levels = 0.05;
end

% sort plist
splist = sort(plist);
levels = levels(:)';
th = zsz(levels);
nlist = numel(splist);
if nargin < 4 || ...
   ~isa(fdrfac, 'double') || ...
    numel(fdrfac) ~= 1 || ...
    isinf(fdrfac) || ...
    isnan(fdrfac) || ...
    fdrfac < 1
    fdrfac = nlist;
end

% divide with FDR factor list for fast find access
nsplist = splist(:)' ./ ((1/ nlist):((1-1/ nlist)/(nlist-1)):1);
fsplist = splist(:)' ./ ((1/fdrfac):((1-1/fdrfac)/(nlist-1)):1);
lastcross = true;
if nargin > 2 && ...
    isa(lnmeth, 'double') && ...
    numel(lnmeth) == 1 && ...
    any((0:3) == lnmeth)
    if lnmeth < 2
        lastcross = false;
    end
    if lnmeth == 1 || ...
        lnmeth == 3
        lnmeth = true;
        th = [th(:), th(:)];
    else
        lnmeth = false;
        th = th(:);
    end
else
    lnmeth = true;
end

% find thresholds
for thc = 1:numel(levels)
    if lastcross
        thi = find(nsplist <= levels(thc));
        if isempty(thi)
            th(thc) = splist(1) / 2.001;
        else
            th(thc) = splist(thi(end));
        end
    else
        thi = find(nsplist > levels(thc));
        if isempty(thi)
            th(thc) = splist(end);
        elseif thi(1) > 1
            th(thc) = splist(thi(1) - 1);
        else
            th(thc) = splist(1) / 2.001;
        end
    end
    if lnmeth
        if lastcross
            thi = find(fsplist <= (levels(thc) / (log(fdrfac) + 0.5772)));
            if isempty(thi)
                th(thc, 2) = splist(1) / 2.001;
            else
                th(thc, 2) = splist(thi(end));
            end
        else
            thi = find(fsplist > (levels(thc) / (log(fdrfac) + 0.5772)));
            if isempty(thi)
                th(thc, 2) = splist(end);
            elseif thi(1) > 1
                th(thc, 2) = splist(thi(1) - 1);
            else
                th(thc, 2) = splist(1) / 2.001;
            end
        end
    end
end
