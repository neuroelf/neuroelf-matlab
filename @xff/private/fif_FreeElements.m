function xo = fif_FreeElements(xo, espec)
% FIF::FreeElements  - memory free elements for file
%
% FORMAT:       fif.FreeElements(espec);
%
% Input fields:
%
%       espec       Element specification, either 1xN double (numbers) or
%                   1x1 struct with one of the fields
%        .block     1x1 double blocktype
%        .range     1x2 double with range of fields to read
%        .type      1x1 double giving the field type
%
% No output fields.
%
% Note: to free all elements, set espec to [Inf].
%
% Using: fifio.

% Version:  v1.1
% Build:    16020310
% Date:     Feb-03 2016, 10:49 AM EST
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

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'fif') || ...
   (~isa(espec, 'double') && ~isstruct(espec)) || isempty(espec)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
espec = espec(:)';

% shortcut to structure
bc = xo.C;
fst = bc.FIFStructure;
nel = numel(fst.Directory);

% second argument is double
if ~isstruct(espec)
    if isinf(espec(1))
        espec = 1:nel;
    else
        if any(espec < 1 | espec > nel | espec ~= fix(espec) | isnan(espec) | isinf(espec))
            error('neuroelf:xff:badArgument', ...
                'Invalid numeric espec argument in call to %s.', mfilename);
        end
    end
else

    % take only first struct
    espec = espec(1);

    % block specification
    if isfield(espec, 'block') && isa(espec.block, 'double') && numel(espec.block) == 1
        bup = fst.BlockLookup;
        bspec = find(bup(1, :) == espec.block);
        espec = zeros(1, 0);
        for buc = bspec
            espec = [espec, bup(2, buc):bup(3, buc)];
        end

    % range specification
    elseif isfield(espec, 'range') && isa(espec.range, 'double') && numel(espec.range) == 2 && ...
        all(espec.range == fix(espec.range) & espec.range > 0 & espec.range <= nel)
        espec = espec.range(1):espec.range(2);

    % type specification
    elseif isfield(espec, 'type') && isa(espec.type, 'double') && numel(espec.type) == 1
        espec = find(fst.Lookup == espec.type);

    % invalid ...
    else
        error('neuroelf:xff:badArgument', 'Invalid struct or bad field contents.');
    end
end

% read elements
fst = ne_methods.fifio(fst, 'freeelem', espec);

% put FIFStructure back into object
bc.FIFStructure = fst;
xo.C = bc;
