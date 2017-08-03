% FUNCTION ne_cm_updatertv: update GLM RunTimeVars
function ne_cm_updatertv(varargin)

% Version:  v1.1
% Build:    16022911
% Date:     Feb-29 2016, 11:57 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2016, Jochen Weber
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
cc = ne_gcfg.fcfg.CM;

% last GLM?
if nargin > 2 && ...
    numel(varargin{3}) == 1 && ...
    isxff(varargin{3}, 'glm')
    glm = varargin{3};
else
    glm = cc.glm;
end

% get options
glm.RunTimeVars.Contrasts = cc.cons;
try
    glm.RunTimeVars.CovariatesData = cat(2, cc.covs{:, 2});
    glm.RunTimeVars.CovariatesNames = cc.covs(:, 1);
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    glm.RunTimeVars.CovariatesData = zeros(glm.NrOfSubjects, 0);
    glm.RunTimeVars.CovariatesNames = cell(0, 1);
end
glm.RunTimeVars.Groups = cc.groups;
if ~cc.usegroups
    glm.RunTimeVars.SubSel = cc.subsel;
    if size(glm.RunTimeVars.SubSels, 1) < 2
        glm.RunTimeVars.SubSels = {'default', cc.subsel};
    else
        glm.RunTimeVars.SubSels{findfirst( ...
            strcmpi(glm.RunTimeVars.SubSels(:, 1), ...
            glm.RunTimeVars.SubSelsSel)), 2} = cc.subsel;
    end
end
glm.RunTimeVars.UseGroups = cc.usegroups;
ne_cm_updateuis(0, 0, cc.glm);

% try to save the RunTimeVars if requested
if nargin > 3 && numel(varargin{4}) == 1 && ...
   (isa(varargin{4}, 'double') || islogical(varargin{4})) && ...
   (true && varargin{4})
    try
        glm.SaveRunTimeVars;
        if ne_gcfg.c.echo
            ne_echo('glm', 'SaveRunTimeVars');
        end
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
end
