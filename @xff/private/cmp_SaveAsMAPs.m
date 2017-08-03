function xo = cmp_SaveAsMAPs(xo, pattern)
% CMP::SaveAsMAPs  - save FMR based component Maps as MAPs
%
% FORMAT:       cmp.SaveAsMAPs(pattern)
%
% Input fields:
%
%       pattern     filename pattern (default: '%s_IC%03d.map')
%                   pattern must contain '%(\d+)?d' and end in '.map'
%                   optionally pattern may contain zero or one '%s'
%
% No output fields.

% Version:  v1.1
% Build:    16020110
% Date:     Feb-01 2016, 10:52 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2011, 2014, 2016, Jochen Weber
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
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'cmp')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
if bc.DocumentType > 1
    error('neuroelf:xff:notSupported', ...
        'Saving MAPs only supported for fmrICA and vtcICA files.');
end
if isempty(bc.Map) || isempty(bc.Map(1).CMPData)
    error('neuroelf:xff:badArgument', 'First component Map must not have empty CMPData.');
end
if nargin < 2 || ~ischar(pattern) || isempty(pattern) || sum(pattern(:) == '%') > 2 || ...
    isempty(regexpi(pattern(:)', '^.*\%\d*d.*\.map$')) || ...
   (sum(pattern(:) == '%') > 1 && isempty(regexpi(pattern(:)', '\%s')))
    pattern = '%s_IC%03d.map';
end
if ~isempty(regexpi(pattern, '\%s'))
    [fp, fn] = fileparts(xo.F);
    fn = {fn};
else
    fn = {};
end

% create required map
map = xff('new:map');
mc = map.C;

% make general settings
cp = bc.Map(1).CMPData;
cs = size(cp);
mc.Type = 0;
mc.DimX = cs(1);
mc.DimY = cs(2);
mc.NrOfSlices = cs(3);
mc.Map(1).Number = 1;
mc.Map(1).Data = single(zeros(cs(1:2)));
mc.Map = mc.Map(ones(1, cs(3)));

% iterate over components
for cc = 1:numel(bc.Map)

    % get CMPData
    cp = bc.Map(cc).CMPData;
    if ~isequal(size(cp), cs)
        continue;
    end

    % set in Map of MAP
    for slc = 1:cs(3)
        mc.Map(slc).Data(:, :) = cp(:, :, slc);
    end

    % try to save as
    map.C = mc;
    newfile = sprintf(pattern, fn{:}, cc);
    try
        aft_SaveAs(map, newfile);
    catch xfferror
        warning('neuroelf:xff:saveError', 'Couldn''t save MAP file ''%s'': %s.', ...
            newfile, xfferror.message);
    end
end

% don't keep object
delete(map);
