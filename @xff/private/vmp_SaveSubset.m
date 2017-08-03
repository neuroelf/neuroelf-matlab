function xo = vmp_SaveSubset(xo, sset, nfilename)
% VMP::SaveSubset  - save a subset of maps to a new file
%
% FORMAT:       vmp.SaveSubset(sset, newfilename);
%
% Input fields:
%
%       sset        1xN list of Maps to save
%       newfilename name for new VMP file
%
% No output fields.

% Version:  v1.1
% Build:    16021320
% Date:     Feb-13 2016, 8:25 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2011, 2014, 2016, Jochen Weber
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
if nargin ~= 3 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'vmp') || ...
   ~isa(sset, 'double') || isempty(sset) || ...
    any(isinf(sset(:)') | isnan(sset(:)') | fix(sset(:)') ~= sset(:)' | sset(:) < 1) || ...
   ~ischar(nfilename) || isempty(nfilename)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;

% correctly get arguments
sset = unique(min(sset(:), bc.NrOfMaps));
if isempty(sset)
    error('neuroelf:xff:badArgument', 'No maps to save remain.');
end

% copy object but only wanted maps
xo2 = aft_CopyObject(xo);
bc2 = xo2.C;
bc2.Map = bc.Map(sset(:)');
bc2.NrOfMaps = numel(bc2.Map);
xo2.C = bc2;

% try to save
try
    aft_SaveAs(xo2, nfilename(:)');
catch xfferror
    warning('neuroelf:xff:internalError', ...
        'Error saving maps to new file: ''%s''.', xfferror.message);
end

% remove object from memory
delete(xo2);
