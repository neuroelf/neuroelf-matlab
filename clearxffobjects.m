function clearxffobjects(olist)
%CLEARXFFOBJECTS  Issues ClearObject call on several objects.
%   CLEARXFFOBJECTS(OLIST) performs the OBJECT.ClearObject call on all
%   valid @xff objects in OLIST. If OLIST is not a cell array, the call
%   does nothing. If any object is the ROOT object, this doesn't lead to
%   an error message.
%
%   See also @XFF.

% Version:  v1.1
% Build:    16012317
% Date:     Jan-23 2016, 5:07 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2016, Jochen Weber
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
if nargin < 1 || ...
   (~iscell(olist) && ~ischar(olist) && ~isa(olist, 'xff')) || isempty(olist)
    return;
end

% for cell array iterate over cells
if iscell(olist)

    % prepare double list for alternative calling syntax
    nlist = char(zeros(0, 24));
    for c = 1:numel(olist)
        if isxff(olist{c})
            delete(olist{c});
        elseif iscell(olist{c})
            clearxffobjects(olist{c});
        elseif ischar(olist{c}) && ndims(olist{c}) < 3 && size(olist{c}, 2) == 24
            nlist = cat(1, nlist, olist{c});
        end
    end

    % clear objects
    if ~isempty(nlist)
        xff(0, 'clearobj', nlist);
    end

% for xff objects
else
    clearxffobjects({olist});
end
