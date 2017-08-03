%XSUBSREF  Non-inclass overloading of SUBSREF.
%   V = XSUBSREF(A, S) retrieves the value V at the sub-locator
%   specified in S.
%
%   While the function behaves like the regular SUBSREF, it is not
%   recommended to call it manually, but it may be used by any user
%   defined class that wishes to access objects of the same class
%   in one of its subfields (which otherwise would break with the
%   dot-syntax):
%
%   Example:
%
%   x = xstruct('field', 'value');
%   y = xstruct('x', x);
%   v = y.x.field
%
%   v =
%
%   value
%
%   Without the implementation *outside* the XSTRUCT class, this
%   would lead to an error message, and the passing on of S to
%   a regular call to subsref *within* the class would try to locate
%   the (hidden) property "field" within the object x, and fail.
%
%   Internally, it simply uses feval with a handle to the general
%   function SUBSASGN, such that this allows class methods to safely
%   overload.
%
%   See also XSTRUCT, SUBSASGN, SUBSREF.

% Version:  v1.0
% Build:    16011717
% Date:     Jan-17 2016, 5:46 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2016, Jochen Weber
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

function [varargout] = xsubsref(a, s)
[varargout{1:nargout}] = subsref(a, s);
