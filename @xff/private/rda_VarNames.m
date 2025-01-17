function [n, t] = rda_VarNames(xo)
% RDA::VarNames  - return a cell array (Nx1) of variable names
%
% FORMAT:       [vnames, vtypes] = obj.MapNames;
%
% No input fields.
%
% Output fields:
%
%       vnames      Nx1 list with variable names
%       vtypes      Nx1 list with variable types

% Version:  v1.1
% Build:    16060121
% Date:     Jun-01 2016, 9:23 PM EST
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

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'rda')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end

% get variables
v = xo.C.Vars;
n = {v.Name};
n = n(:);

% get types as well
if nargout > 1
    t = {v.Data};
    t = t(:);
    for vc = 1:numel(t)
        if isstruct(t{vc}) && numel(t{vc}) == 1 && isfield(t{vc}, 'Type')
            t{vc} = t{vc}.Type;
        else
            t{vc} = class(t{vc});
        end
    end
end
