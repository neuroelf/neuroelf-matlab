function [p, pcol] = glm_SubjectPredictors(xo)
% GLM::SubjectPredictors  - return list of subject predictor names
%
% FORMAT:       [spreds, spredcol] = glm.SubjectPredictors;
%
% No input fields.
%
% Output fields:
%
%       spreds      subject predictor names list (Px1)
%       spredcol    subject predictor colors (Px3 RGB code)
%
% Using: findfirst.

% Version:  v1.1
% Build:    16020409
% Date:     Feb-04 2016, 9:37 AM EST
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

% neuroelf library
global ne_methods;

% check arguments
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'glm')
    error('neuroelf:xff:badArgument', 'Invalid object handle in call.');
end

% get content
bc = xo.C;

% get ALL subjects' predictor names
p = {bc.Predictor(:).Name2};

% allow for empty returns (on new objects)
if isempty(p)
    pcol = zeros(0, 3);
    return;
end

% confound name (for fixed name, without subject)
pl = p{end};

% which to keep
if bc.ProjectTypeRFX > 0 || bc.SeparatePredictors == 2

    % get first subject's ID
    fs = p{1}(1:ne_methods.findfirst(p{1} == ':'));

    % replace that with nothing
    p = strrep(p, [fs ' '], '');

    % get to keep
    keep = cellfun('isempty', regexp(p, '^Subject\s+.*:\s+')) & ~strcmpi(p, pl);

    % FFX ?
    if bc.ProjectTypeRFX < 1
        keep = keep & cellfun('isempty', regexp(p, '^Study'));
    end

% regular FFX
else

    % study based FFX
    if bc.SeparatePredictors == 1
        keep = ~cellfun('isempty', regexp(p, '^Study 1:'));
        keepco = 1 + ne_methods.findfirst(diff(double(keep)));
        keep(keepco:end) = false;
        p(keep) = regexprep(p(keep), '^Study\s+1:\s+', '');
    else
        keep = cellfun('isempty', regexp(p, '^Study'));
    end
end

% keep only those
p = p(keep);
p = p(:);

% colors also?
if nargout > 1
    pcol = {bc.Predictor(:).RGB};
    pcol = cat(3, pcol{:});
    pcol = permute(pcol(1, :, keep), [3, 2, 1]);
end
