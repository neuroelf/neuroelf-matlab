% FUNCTION ne_glmgeneratemdm: generate MDM file
function [varargout] = ne_glmgeneratemdm(varargin)

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:33 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, 2016, Jochen Weber
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

% global variable
global ne_gcfg;
cc = ne_gcfg.fcfg;

% preset output
varargout = cell(1, nargout);

% only valid for certain combinations
if numel(cc.StatsVar) ~= 1 || ...
   ~isxff(cc.StatsVar, 'glm') || ...
    cc.StatsVar.ProjectTypeRFX ~= 1
    return;
end

% try call
mdm = {[]};
try
    mdm{1} = cc.StatsVar.GenerateMDM;
    mdm{1}.SaveAs;
catch ne_eo;
    clearxffobjects(mdm);
    uiwait(warndlg(sprintf( ...
        'An error occurred while trying to create MDM from GLM:\n%s', ...
        ne_eo.message), 'modal'));
    return;
end

% output
if nargout > 0
    varargout{1} = mdm{1};
else
    clearxffobjects(mdm);
end
