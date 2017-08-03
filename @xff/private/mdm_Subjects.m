function s = mdm_Subjects(xo, full)
% MDM::Subjects  - return list of subjects of multi-subject MDM
%
% FORMAT:       subjects = mdm.Subjects([full]);
%
% Input fields:
%
%       full        flag, if true, do not return unique list
%
% Output fields:
%
%       subjects    subjects list (Sx1 cell array)

% Version:  v1.1
% Build:    16012617
% Date:     Jan-26 2016, 5:56 PM EST
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

% check arguments
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'mdm')
    error('neuroelf:xff:badArgument', 'Invalid object handle in call.');
end
if nargin < 2 || ~islogical(full) || numel(full) ~= 1
    full = false;
end
bc = xo.C;
s = bc.XTC_RTC(:, 1);
for sc = 1:numel(s)
    [null, s{sc}] = fileparts(s{sc});
end
s = regexprep(s, '[_\.].*$', '');
if full
    return;
end
[su, sui] = unique(s);
s = s(sort(sui));
