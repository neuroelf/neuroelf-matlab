function xo = aft_ClearObject(xo)
% AFT::ClearObject  - remove object from global storage
%
% FORMAT:       obj.ClearObject;
%
% No input / output fields.
%
% TYPES: ALL

% Version:  v1.1
% Build:    16022911
% Date:     Feb-29 2016, 11:29 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2016, Jochen Weber
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
if nargin ~= 1 || numel(xo) ~= 1 || ~xffisobject(xo, true) || xo.L(1) == 'X'
    error('neuroelf:xff:badArgument', 'Invalid xff object in call.');
end

% save run time vars
rtv = xo.C.RunTimeVars;
if isfield(rtv, 'AutoSave') && islogical(rtv.AutoSave) && numel(rtv.AutoSave) == 1 && rtv.AutoSave && ...
    isfield(xo.H, 'RunTimeVarsSaved') && islogical(xo.H.RunTimeVarsSaved) && ...
    numel(xo.H.RunTimeVarsSaved) == 1 && ~xo.H.RunTimeVarsSaved && ~isempty(xo.F)
    try
        aft_SaveRunTimeVars(xo);
    catch xfferror
        warning('neuroelf:xff:SaveRunTimeVarsOnClear', xfferror.message);
    end
end

% call clear BEFORE delete (so it can still function)
xffclear(xo.L);

% call delete
delete(xo);
