function rstruct = renamefields(istruct, ffrom, fto)
% renamefields  - renames fields in a struct
%
% FORMAT:       renamed = renamefields(istruct, ffrom, fto)
%
% Input fields:
%
%       istruct     input structure array
%       ffrom       from field names cell array
%       fto         to field names cell array
%
% See also rmfield.

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
if nargin < 3 || ...
   ~isstruct(istruct) || ...
   ~iscell(ffrom) || ...
   ~iscell(fto) || ...
    isempty(istruct) || ...
    numel(ffrom) ~= numel(fto)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or too few arguments passed.' ...
    );
end

% get input sizes
flds = fieldnames(istruct);
nfld = length(flds);
sstr = size(istruct);
nstr = prod(sstr);

% get field sizes
ffrom = ffrom(:);
fto   = fto(:);
nflds = length(ffrom);

% remove unknown/bad fields
for fc = nflds:-1:1
    if ~ischar(ffrom{fc}) || ...
       ~ischar(fto{fc})
        ffrom(fc) = [];
        fto(fc)   = [];
        continue;
    end
    ffrom{fc} = ffrom{fc}(:)';
    fto{fc}   = fto{fc}(:)';
    if ~isrealvarname(ffrom{fc}) || ...
       ~isrealvarname(fto{fc}) || ...
        strcmp(ffrom{fc}, fto{fc}) || ...
       ~isfield(istruct, ffrom{fc})
        ffrom(fc) = [];
        fto(fc)   = [];
        continue;
    end
end
nflds = length(ffrom);
if nflds == 0
    rstruct = istruct;
    return;
end

% build lookup structs
fstruct = struct;
for fc = 1:nflds
    fstruct.(ffrom{fc}) = fc;
end
nstruct = struct;
for fc = 1:nfld
    if isfield(fstruct, flds{fc})
        nstruct.(flds{fc}) = true;
    else
        nstruct.(flds{fc}) = false;
    end
end

% build output structure
rstruct = struct;

% build output size argument if needed
if any(sstr > 1)
    osize = sprintf('%.0f,', sstr);
    eval(['rstruct(' osize(1:end-1) ')=rstruct;']);
end

% do the work
for sc = 1:nstr
    for fc = 1:nfld
        if nstruct.(flds{fc})
            rstruct(sc).(fto{fstruct.(flds{fc})}) = istruct(sc).(flds{fc});
        else
            rstruct(sc).(flds{fc}) = istruct(sc).(flds{fc});
        end
    end
end
