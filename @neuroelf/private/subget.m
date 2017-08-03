function pstruct = subget(hnd,props)
% subget  - get a list of properties as a struct
%
% this function simplifies getting multiple properties from
% MATLAB's graphics handle objects.
%
% FORMAT:       propstruct = subget(handle, props)
%
% Input fields:
%
%       handle      Nx1 graphic handles
%       props       either a single string or a cell array or
%                   strings naming properties of the handles
%
% Output fields:
%
%       propstruct  Nx1 structure array with property contents
%
% See also get.

% Version:  v1.0
% Build:    15040310
% Date:     Apr-03 2015, 10:11 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2015, Jochen Weber
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

% enough arguments ?
if nargin < 2
    error( ...
        'neuroelf:TooFewArguments',...
        'Too few arguments. Try ''help %s''.',...
        mfilename ...
    );
end

% not all inputs are valid handles
if isempty(hnd) || ...
   ~all(ishandle(hnd))
    error( ...
        'neuroelf:BadArgument',...
        'Invalid handle argument.' ...
    );
end

% cell-ify char/string props if necessary
if ischar(props)
    props = {props(:)'};
end

% not a valid cell?
if ~iscell(props)
    error( ...
        'neuroelf:BadArgument',...
        'Invalid props argument.' ...
    );
end

% prepare for getting values
hnd     = hnd(:);
nhnd    = size(hnd, 1);
pstruct = cell2struct(cell(1, 1, 0), {}, 3);

% iterate over properties
for cc = 1:numel(props)
    
    % only keep valid ones (tested on first handle only!)
    try
        get(hnd(1), props{cc});
        pstruct.(props{cc}) = [];
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
    end
end

% use valid properties
props = fieldnames(pstruct);
if isempty(props)
    error( ...
        'neuroelf:NoValidProperty',...
        'No valid property names given.' ...
    );
end

% prepare
props = props(:)';
pstruct(1) = [];

% and get values
values = get(hnd,props);
for cc = 1:size(props,2)
    [pstruct(1:nhnd, 1).(props{cc})] = deal(values{:, cc});
end
