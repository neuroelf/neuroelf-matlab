function xo = smp_Combine(xo, option, ospec)
% SMP::Combine  - combines maps within an SMP object
%
% FORMAT:       [smp = ] smp.Combine(option, ospec);
%
% Input fields:
%
%       option      string (what to do), either of
%                   - average:        calculate average
%       ospec       1x1 struct with appropriate fields
%       .group      selection of maps to apply operation on
%       .name       name of combined map (if any)
%
% Output fields:
%
%       smp         object with altered maps or one added, combined map
%
% Note:  map options (thresholds, etc.) are copied from first map in
%        ospec.group

% BEGIN to do:
%                   - conj, max, min: conjunction, minimum or maximum map
%                   - correl:         2-group correlation (size match!)
%                   - func:           apply func and store result
%                   - onesample:      one-sample t-test (group > 0)
%                   - paired:         paired t-test (size match!)
%                   - twosample:      two-sample t-test (group > group2)
%       .func       if given, perform this function on each map before
%                   combination
%                   thereby '@' is replaced by the map values
%                   - '@(@ < 0) = 0' would set all values < 0 to 0.
%                   - '@ = @ - mean(@)' would remove the mean from the map
%                   and '$x' is replaced by the map values of map x
%                   - '@ = @ - $1' subtract map one from current map
%       .func2      function applied on maps of second selection
%       .funcname   pattern to use to alter map names (must contain %s)
%       .funcname2  pattern of second selection
% NOTE:             functions can be daisy-chained by ';'
%       .group2     second selection of maps (for binary operations)
%       .mask       1x2 double array for masking (default: [-Inf, Inf])
%       .mask2      masking option of second selection
% END to do:

% Version:  v1.1
% Build:    16021110
% Date:     Feb-11 2016, 10:21 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2014, 2016, Jochen Weber
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
if nargin > 2 && isa(ospec, 'double')
    ospec = struct('group', ospec);
end
if nargin < 3 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'smp') || ...
   ~ischar(option) || ~any(strcmpi(option(:)', {'average', 'avg', 'conj', ...
        'conjunction', 'cor', 'correl', 'func', 'function', 'max', 'maximum', ...
        'min', 'minimum', 'onesample', 'paired', 'twosample'})) || ...
   ~isstruct(ospec) || numel(ospec) ~= 1 || ~isfield(ospec, 'group') || ...
   ~isa(ospec.group, 'double') || isempty(ospec.group) || ...
    any(isinf(ospec.group(:)) | isnan(ospec.group(:)) | ...
        ospec.group(:) ~= fix(ospec.group(:)) | ospec.group(:) < 1)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
if any(ospec.group(:) > numel(bc.Map))
    error('neuroelf:xff:badArgument', 'Specified map(s) out of bounds.');
end
option = lower(option(:)');
ospec.group = ospec.group(:)';
if ~isfield(ospec, 'group2') || ~isa(ospec.group2, 'double') || isempty(ospec.group2) || ...
    any(isinf(ospec.group2(:)) | isnan(ospec.group2(:)) | ...
        ospec.group2(:) ~= fix(ospec.group2(:)) | ospec.group2(:) < 1 | ...
        ospec.group2(:) > numel(bc.Map))
    ospec.group2 = [];
else
    ospec.group2 = ospec.group2(:)';
end
if ~isfield(ospec, 'name') || ~ischar(ospec.name) || isempty(ospec.name)
    ospec.name = '';
else
    ospec.name = ospec.name(:)';
end

% check integrity
gsel1 = ospec.group;
gsel2 = ospec.group2;
gsel = [gsel1, gsel2];
mtype = bc.Map(gsel(1)).Type;
numv = numel(bc.Map(gsel(1)).SMPData);
for mc = gsel(2:end)
    if bc.Map(mc).Type ~= mtype
        error('neuroelf:xff:badArgument', 'Invalid combination of map types.');
    end
    if numel(bc.Map(mc).SMPData) ~= numv
        error('neuroelf:xff:badObject', 'Map sizes mismatch. Invalid object!');
    end
end

% what to do
switch (option)

    % average map of group 1
    case {'average', 'avg'}

        % get first map
        mvals = bc.Map(gsel1(1)).SMPData(:);
        nummp = numel(gsel1);

        % for all but CC maps
        if mtype ~= 3

            % add 2nd to n-th map
            for mc = gsel1(2:end)
                mvals = mvals + bc.Map(mc).SMPData(:);
            end

            % divide by number of maps
            mvals = mvals / nummp;

        % CC maps
        else

            % split map values
            rvals = mvals - floor(mvals);
            lvals = round(mvals - rvals) / 1000;

            % procede with 2nd to n-th map
            for mc = gsel1(2:end)

                % get values
                avals = bc.Map(mc).SMPData(:);

                % split r-values
                arvals = avals - floor(avals);

                % add to 1st map values
                rvals = rvals + arvals;
                lvals = lvals + round(avals - arvals) / 1000;
            end

            % build avarage
            mvals = min(1-eps, max(0, rvals / nummp)) + 1000 * round(lvals / nummp);
        end

        % add map
        bc.Map(end+1) = bc.Map(gsel1(1));
        if isempty(ospec.name)
            bc.Map(end).Name = ['Average map:' sprintf(' %d', gsel1)];
        else
            bc.Map(end).Name = ospec.name;
        end
        bc.Map(end).SMPData = mvals;
        bc.NrOfMaps = numel(bc.Map);

%    case {'conj', 'conjunction'}
%    case {'max', 'maximum'}
%    case {'min', 'minimum'}
%    case {'cor', 'correl'}
%    case {'func', 'function'}
%    case {'onesample'}
%    case {'paired'}
%    case {'twosample'}
    otherwise
end

% set content back
xo.C = bc;
