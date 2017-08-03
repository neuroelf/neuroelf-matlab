function pv = mappvalue(v, map, iv, ot)
% mappvalue  - convert stats to pvalue (signed)
%
% FORMAT:       pv = mappvalue(v, map [, iv [, ot]])
%
% Input fields:
%
%       v           value (e.g. from Map.VMPData or Map.SMPData)
%       map         1x1 struct with minimum fields
%        .DF1       d.f.1 parameter
%        .DF2       d.f.2 parameter (for F-stats)
%        .Type      numeric type, one of
%                   1 (t-stat)
%                   2/3 (r-stat)
%                   4 (F-stat)
%                   12 (z-stat)
%       iv          inverse operation flag (default: false)
%       ot          one-tailed flag (default: false)
%
% Output fields:
%
%       pv          p-value
%
% Note: by default, the p-values will be two-tailed and signed to
%       allow easy use of max/min computations; to get the

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
if nargin < 2 || ...
   (~isa(v, 'single') && ...
    ~isa(v, 'double')) || ...
    numel(map) ~= 1 || ...
   ~isstruct(map) || ...
   ~isfield(map, 'DF1') || ...
   ~isa(map.DF1, 'double') || ...
    numel(map.DF1) ~= 1 || ...
    isinf(map.DF1) || ...
    isnan(map.DF1) || ...
    map.DF1 < 1 || ...
   ~isfield(map, 'Type') || ...
   ~isa(map.Type, 'double') || ...
    numel(map.Type) ~= 1 || ...
    isinf(map.Type) || ...
    isnan(map.Type) || ...
   ~any([1:4, 12] == map.Type) || ...
   (map.Type == 4 && ...
    (~isfield(map, 'DF2') || ...
     ~isa(map.DF2, 'double') || ...
      numel(map.DF2) ~= 1 || ...
      isinf(map.DF2) || ...
      isnan(map.DF2) || ...
      map.DF2 < 1))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end

% inverse flag
if nargin < 3 || ...
   ~islogical(iv) || ...
    numel(iv) ~= 1
    iv = false;
end
if nargin < 4 || ...
   ~islogical(ot) || ...
    numel(ot) ~= 1
    ot = false;
end

% don't allow negative values for one-tailed stats or F-stats
if (ot && ...
    iv) || ...
    map.Type == 4
    if ~iv
        v(v < 0) = 0;
    else
        v(v <= 0) = 1;
    end
end

% for z-stats, simply set DF1 to 1e7
if map.Type == 12
    map.DF1 = 1e7;
end

% regular mode
if ~iv

    % depending on stats type
    switch (map.Type)

        % t-stats
        case {1, 12}
            pv = 2 * (sign(v) .* sdist('tcdf', -abs(double(v)), map.DF1));

        % r-stats
        case {2, 3}

            % remove CC-lag portion
            if map.Type == 3
               v = v - floor(v);
            end
            pv = sign(v) .* correlpvalue(double(v), map.DF1 + 2);

        % F-stats
        case {4}
            pv = (v >= 0) .* sdist('fcdf', double(v), map.DF1, map.DF2, true);
    end

    % one tailed
    if ot && ...
        map.Type ~= 4
        ps = (pv >= 0);
        pv(ps) = 1 - 0.5 * pv(ps);
        pv(~ps) = -0.5 * pv(~ps);
    end

% inverse mode
else

    % from one-tailed stats
    if ot && ...
        map.Type ~= 4
        ps = (v >= 0.5);
        v(ps) = 2 - 2 * v(ps);
        v(~ps) = -2 * v(~ps);
    end

    % depending on stats type
    switch (map.Type)

        % t-stats
        case {1, 12}
            pv = -sign(v) .* sdist('tinv', 0.5 * abs(double(v)), map.DF1);

        % r-stats
        case {2, 3}
            pv = correlinvtstat( ...
                -sign(v) .* sdist('tinv', 0.5 * abs(double(v)), map.DF1), map.DF1 + 2);

        % F-stats
        case {4}
            pv = (v >= 0) .* sdist('finv', double(v), map.DF1, map.DF2, true);
    end
end
