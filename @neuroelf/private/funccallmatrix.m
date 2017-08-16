function fcm = funccallmatrix(mfiles)
%FUNCCALLMATRIX  Create a boolean matrix for function calls in files.
%   FCM = FUNCCALLMATRIX(MFILES) reads each file in MFILES and records
%   all occurrences of the other files. If the (local) names in MFILES
%   are not unique, all potential matches will be selected. The returned
%   FCM variable is a FxF boolean matrix such that if MFILES{FROM}
%   contains a call to MFILES{to} FCM(FROM, TO) will be true.
%
%   Example:
%
%   mfiles = findfiles(toolboxdir, '*.m');
%   fcm = funccallmatrix(mfiles);

% Version:  v1.1
% Build:    17081521
% Date:     Aug-15 2017, 9:07 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2017, Jochen Weber
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
if nargin < 1 || ~iscell(mfiles) || isempty(mfiles) || ...
   ~all(cellfun(@ischar, mfiles(:))) || any(cellfun('isempty', mfiles(:)))
    error('neuroelf:general:badOrMissingInput', 'Bad or missing input.');
end

% load files
try
    nfiles = numel(mfiles);
    mfiles = mfiles(:);
    mfcont = cell(nfiles, 1);
    for fc = 1:nfiles
        mfcont{fc} = asciiread(mfiles{fc});
    end
    fcm = false(nfiles, nfiles);
catch ne_eo;
    rethrow(ne_eo);
end

% get short names
[mfolders, mfilesh] = mfileparts(mfiles);

% for now, use brute force
for fc = 1:nfiles
    for sc = 1:nfiles
        if sc == fc
            continue;
        end

        % look for the symbol in the file
        symfound = regexpi(mfcont{fc}, ['[^a-zA-Z_]' mfilesh{sc} '\s*\(']);
        if isempty(symfound)
            symfound = regexpi(mfcont{fc}, ['\@' mfilesh{sc} '[^a-zA-Z_0-9]']);
        end
        if ~isempty(symfound)
            symassigned = regexpi(mfcont{fc}, ['[^a-zA-Z_]' mfilesh{sc} '\s*=']);
            if isempty(symassigned) || symassigned(1) > symfound(1)
                fcm(fc, sc) = true;
            end
        end
    end
end