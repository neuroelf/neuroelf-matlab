function s = glm_Subjects(xo, full)
% GLM::Subjects  - return list of subjects of multi-subject GLM
%
% FORMAT:       subjects = glm.Subjects([full]);
%
% Input fields:
%
%       full        flag, if true, returns subject IDs for each study
%
% Output fields:
%
%       subjects    subjects list (Sx1 cell array)

% Version:  v1.1
% Build:    16020316
% Date:     Feb-03 2016, 4:00 PM EST
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
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'glm')
    error('neuroelf:xff:badArgument', 'Invalid object handle in call.');
end
if nargin < 2 || ~islogical(full) || numel(full) ~= 1
    full = false;
end
bc = xo.C;
if isfield(xo.H, 'SubjectIDs') && iscell(xo.H.SubjectIDs) && numel(xo.H.SubjectIDs) == numel(bc.Study)
    s = xo.H.SubjectIDs;
else
    s = {bc.Study(:).NameOfAnalyzedFile};
    s = s(:);
    for sc = 1:numel(s)
        [p, s{sc}] = fileparts(s{sc});
        s{sc} = regexprep(s{sc}, '^([^_]+)_.*$', '$1');
    end
    xo.H.SubjectIDs = s;
end
if full
    return;
end
[su, sui] = unique(s);
if bc.ProjectTypeRFX > 0 && numel(su) ~= bc.NrOfSubjects
    warning('neuroelf:xff:internalError', 'NrOfSubjects does not match with unique subject IDs.');
end
s = s(sort(sui));
