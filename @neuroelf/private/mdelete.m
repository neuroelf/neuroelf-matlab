function succ = mdelete(files)
% mdelete  - delete multiple files
%
% FORMAT:       succ = mdelete(files)
%
% Input fields:
%
%       files       cell array with file names (or cell-of-cell arrays)
%
% Output fields:
%
%       succ        success (boolean flag, for cell array)

% Version:  v0.9d
% Build:    14090510
% Date:     Sep-05 2014, 10:51 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2014, Jochen Weber
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

% requires cell array
if nargin < 1 || ...
   ~iscell(files)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end

% prepare output
succ = false(size(files));

% run
for fc = 1:numel(files)
    
    % correct contents (single file)
    if ischar(files{fc}) && ...
       ~isempty(files{fc})
        try
            delete(files{fc}(:)');
            succ(fc) = true;
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
        end

    % possibly cell with files
    elseif iscell(files{fc}) && ...
       ~isempty(files{fc})

        % pass on and use all()
        succ(fc) = all(mdelete(files{fc}));
    end
end
