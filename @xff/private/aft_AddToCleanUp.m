function xo = aft_AddToCleanUp(xo, cleanup)
% AFT::AddToCleanUp  - add a statement to CleanUp handle
%
% FORMAT:       [obj = ] obj.AddToCleanUp(statement)
%
% Input fields:
%
%       statment    1xN char with a valid statement (eval'ed in BASE)
%
% Output fields:
%
%       obj         reference to input object
%
% TYPES: ALL
%
% Using: checksyntax.

% Version:  v1.1
% Build:    16012420
% Date:     Jan-24 2016, 8:10 PM EST
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

% only valid for single file
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true) || xo.L(1) == 'X' || ...
   ~ischar(cleanup) || isempty(cleanup) || ~isempty(ne_methods.checksyntax(cleanup(:)'))
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end

% then add to cleanup
if ~isfield(xo.H, 'CleanUp')
    xo.H.CleanUp = {cleanup(:)'};
else
    xo.H.CleanUp{end+1} = cleanup(:)';
end
