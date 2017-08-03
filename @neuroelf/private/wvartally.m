function [t, varargout] = wvartally(t, s, w, w2)
%WVARTALLY  Tally (weighted) scores and compute mean (difference).
%   T = WVARTALLY(T, S) tallies score S into total T. To begin a tally
%   operation, pass in an empty array as T. Subsequent calls then should
%   use the T as it is returned by WVARTALLY. S can be any size, as long as
%   it is not empty, but must be numeric.
%
%   T = WVARTALLY(T, S, W) tallies score S using weight(s) W, with W being
%   either 1x1 or the same size as S and of either single or double class.
%
%   T = WVARTALLY(T, S, -1, W) is using multiple negative weights.
%
%   [M, V, SE, DF] = WVARTALLY(T) returns the mean, var, and SE, and DF.
%   In the case of a difference, the degrees of freedom will be corrected
%   for difference in variance (according to a two-sample t-test).

% Version:  v1.1
% Build:    16050618
% Date:     May-06 2016, 6:38 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2016, Jochen Weber
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
if nargin < 1 || ~isstruct(t) || numel(t) ~= 1 || ~isfield(t, 'IsTally') || ...
   ~islogical(t.IsTally) || numel(t.IsTally) ~= 1 || ~t.IsTally
    t = struct('IsTally', true, ...
        'NN', 0, 'NS', 0, 'NSS', 0, 'NSW', 0, 'NWM', 0, 'NWS', 0, 'NWW', 0, 'NWX', 0, ...
        'PN', 0, 'PS', 0, 'PSS', 0, 'PSW', 0, 'PWM', 0, 'PWS', 0, 'PWW', 0, 'PWX', 0, ...
        'S',  0, 'SS', 0, 'Size', []);
end

% output check
if nargout > 1
    varargout = cell(1, nargout - 1);
end

% add one tally to another (using variance as weight)
if nargin > 1  && isstruct(s) && numel(s) == 1 && isfield(s, 'IsTally') && ...
    islogical(s.IsTally) && numel(s.IsTally) == 1 && s.IsTally

    % compute mean and variance
    [m, v] = wvartally(s);
    v = 1 ./ v;

    % tally using (1/variance), allowing single weight
    if nargin < 3 || ~isa(w, 'double')
        t = wvartally(t, m, v);
    elseif numel(w) == 1
        t = wvartally(t, m, w, v);
    elseif isequal(size(w), size(v))
        t = wvartally(t, m, w .* v);
    else
        t = wvartally(t, m, v);
    end

    % return
    return;
end

% no more inputs
if nargin < 2 || ~isnumeric(s) || isempty(s)

    % compute mean
    to = t;

    % no negative weights (straight weighted mean)
    if to.NN == 0

        % compute 1 over sum-of-weights (as scaling factor)
        w = 1 ./ to.PSW;

        % compute mean (and store as "t" = first output)
        t = w .* to.PS;

        % variance
        if nargout > 1
            
            % compute bias-corrected factor for variance
            pw = ((to.PN - 1) / to.PN) .* to.PSW;

            % compute and store variance in first varargout
            varargout{1} = (to.PSS - w .* (to.PS .* to.PS)) ./ pw;

            % standard error
            if nargout > 2

                % Cochran W.G. (1977) Sampling Techniques 3rd Ed., see also
                % Gatz, D.F. & Smith, L., "The Standard Error of a Weighted
                % Mean Concentration - Bootstrapping vs. other Methods"
                % (1995) Atmospheric Environment, Vol. 29, No. 11
                varargout{2} = sqrt((to.PN ./ ...
                    ((to.PN - 1) .* to.PSW .* to.PSW)) .* ...
                    ((to.PWW - to.PS .* to.PS ./ to.PN) - ...
                    2 .* t .* (to.PWX - to.PSW .* to.PS ./ to.PN) + ...
                    t .* t .* (to.PWS - to.PSW .* to.PSW ./ to.PN)));

                % DF
                if nargout > 3
                    varargout{3} = to.PSW ./ to.PWM;
                end
            end
        end

    % difference of weighted sums
    else

        % compute 1 over sum-of-weights (as scaling factors)
        pw = 1 ./ to.PSW;
        nw = 1 ./ to.NSW;

        % compute difference of mean (and store as "t" = first output)
        p = (pw .* to.PS);
        n = (nw .* to.NS);
        t = p - n;

        % variance
        if nargout > 1
            
            % compute bias-corrected factor for variance
            pw = (to.PN / (to.PN - 1)) .* pw;
            nw = (to.NN / (to.NN - 1)) .* nw;

            % compute and store variance in first varargout
            varargout{1} = (pw + nw) .* ... 
                ((to.PSS - pw .* (to.PS .* to.PS)) + (to.NSS - nw .* (to.NS .* to.NS)));

            % standard error
            if nargout > 2

                % see above
                vr1 = (to.PN ./ ((to.PN - 1) .* to.PSW .* to.PSW)) .* ...
                    ((to.PWW - to.PS .* to.PS ./ to.PN) - ...
                    2 .* p .* (to.PWX - to.PSW .* to.PS ./ to.PN) + ...
                    p .* p .* (to.PWS - to.PSW .* to.PSW ./ to.PN));
                vr2 = (to.NN ./ ((to.NN - 1) .* to.NSW .* to.NSW)) .* ...
                    ((to.NWW - to.NS .* to.NS ./ to.NN) - ...
                    2 .* n .* (to.NWX - to.NSW .* to.NS ./ to.NN) + ...
                    n .* n .* (to.NWS - to.NSW .* to.NSW ./ to.NN));
                varargout{2} = sqrt(vr1 + vr2);
                if nargout > 3

                    % correct degrees of freedom for difference in variance
                    varargout{3} = ((vr1 + vr2) .^ 2) ./ ...
                        (vr1 .* vr1 ./ (to.PN - 1) + vr2 .* vr2 ./ (to.NN - 1));
                end
            end
        end
    end

    % return
    return;
end

% ensure S is double (we know S is given and valid!)
s = double(s);
sz = size(s);
if ~isempty(t.Size) && ~isequal(sz, t.Size)
    error('neuroelf:general:badArgument', 'S must match in size.');
elseif isempty(t.Size)
    t.Size = sz;
end

% check weights W
if nargin < 3 || isempty(w) || (~isa(w, 'double') && ~isa(w, 'single'))
    w = 1;
elseif numel(w) ~= 1 && ~isequal(size(w), size(s))
    error('neuroelf:general:badArgument', 'W must be 1x1 or match with size(S).');
else
    w = double(w);
end

% add to lowest level
if numel(w) > 1 || w >= 0
    if nargin > 3 && isa(w2, 'double') && isequal(size(w2), size(s))
        if isequal(w, 1)
            w = w2;
        else
            w = w .* w2;
        end
    end
    t.PN = t.PN + 1;
    t.PS = t.PS + w .* s;
    t.PSS = t.PSS + w .* s .* s;
    t.PSW = t.PSW + w;
    t.PWM = max(t.PWM, w);
    t.PWS = t.PWS + w .* w;
    t.PWW = t.PWW + w .* w .* s .* s;
    t.PWX = t.PWX + w .* w .* s;
    t.S = t.S + s;
    t.SS = t.S + s .* s;
else
    if nargin > 3 && isa(w2, 'double') && isequal(size(w2), size(s))
        if isequal(w, -1)
            w = -w2;
        else
            w = w .* w2;
        end
    end
    w = -w;
    t.NS = t.NS + w .* s;
    t.NN = t.NN + 1;
    t.NSS = t.NSS + w .* s .* s;
    t.NSW = t.NSW + w;
    t.NWM = max(t.NWM, w);
    t.NWS = t.NWS + w .* w;
    t.NWW = t.NWW + w .* w .* s .* s;
    t.NWX = t.NWX + w .* w .* s;
end
