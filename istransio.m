function isf = istransio(htio, varargin)
% istransio  - check (and validate) object
%
% FORMAT:       isf = istransio(hfile [, valid])
%
% Input fields:
%
%       hfile       MxN argument check for class
%       valid       if given and true, perform validation
%
% Output fields:
%
%       isf         either plain false or logical array of input
%                   size with check result (for validity check)

% Version:  v0.9d
% Build:    14082210
% Date:     Aug-22 2014, 10:51 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2014, Jochen Weber
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

% base argument check
if nargin < 1
    error( ...
        'neuroelf:TooFewArguments', ...
        'At least one input argument is required.' ...
    );
end
chstrict = false;
if nargin > 1 && ...
    numel(varargin{1}) == 1 && ...
   (isnumeric(varargin{1}) || ...
    islogical(varargin{1}))
    if varargin{1}
        chstrict = true;
    end
end

% class check
if isa(htio, 'transio')
    stio = struct(htio);
    isf = true(size(stio));
else
    isf = false;
    return;
end

% validation ?
if chstrict

    % make struct and check fields
    sflds = fieldnames(stio);
    if length(sflds) ~= 8 || ...
       ~all(strcmp(sflds, ...
           {'FileName'; ...
            'LittleND'; ...
            'DataType'; ...
            'TypeSize'; ...
            'IOOffset'; ...
            'DataDims'; ...
            'ReadOnly'; ...
            'IOBuffer'}))
        isf(:) = false;
        return;
    end
end
