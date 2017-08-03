% FUNCTION ne_cm_copycons: copy contrasts from other GLM
function ne_cm_copycons(varargin)

% Version:  v1.0
% Build:    16011018
% Date:     Jan-10 2016, 6:53 PM EST
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

% global variable
global ne_gcfg;
ch = ne_gcfg.h.CM;

% ask for filename
[glmfile, glmpath] = uigetfile( ...
    {'*.glm',       'General Linear Model files (*.glm)'}, ...
     'Please select the GLM file containing the contrasts...', ...
     'MultiSelect', 'off');
if isequal(glmfile, 0) || ...
    isequal(glmpath, 0) || ...
    isempty(glmfile)
    return;
end
glmfile = [glmpath glmfile];

% try reading contrasts
try
    cfp = ch.CMFig.Pointer;
    ch.CMFig.Pointer = 'watch';
    drawnow;
    glm = xff(glmfile);
    conpred = glm.SubjectPredictors;
    cons = glm.RunTimeVars.Contrasts;
    concols = glm.RunTimeVars.ContrastColors;
    glm.ClearObject;
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    uiwait(warndlg('Invalid GLM file.', 'NeuroElf GUI - Info', 'modal'));
    ch.CMFig.Pointer = cfp;
    return;
end

% no contrasts defined
if isempty(cons)
    uiwait(warndlg('Selected GLM does not contain contrasts.', 'modal'));
    ch.CMFig.Pointer = cfp;
    return;
end

% get config
cc = ne_gcfg.fcfg.CM;
glm = cc.glm;

% predictors cannot be resolved?
glmpred = glm.SubjectPredictors;
glmpred(end) = [];
pmatch = multimatch(glmpred, conpred);
if any(pmatch < 1) || ...
    numel(pmatch) ~= numel(unique(pmatch))
    uiwait(warndlg('Couldn''t resolve all required predictor names.', 'modal'));
    ch.CMFig.Pointer = cfp;
    return;
end

% which contrasts to add
[consel, selok] = listdlg( ...
    'ListString', cons(:, 1), ...
    'SelectionMode', 'multiple', ...
    'ListSize', [min(640, 40 + 10 * size(char(cons(:, 1)), 2)), ...
        min(300, 60 + 10 * size(cons, 1))], ...
    'InitialValue', (1:size(cons, 1))', ...
    'Name', 'Please select the contrasts to add...');
if isempty(consel) || ...
   ~isequal(selok, 1)
    ch.CMFig.Pointer = cfp;
    return;
end

% update RTV first
ne_cm_updatertv;

% compile new weights
conweights = cell(numel(consel), 1);
for conc = 1:numel(consel)
    conweights{conc} = lsqueeze(cons{consel(conc), 2}(pmatch));
end

% then add to GLM
glm.RunTimeVars.Contrasts(end+1:end+numel(consel), :) = ...
    [cons(consel, 1), conweights];
glm.RunTimeVars.ContrastColors(end+1:end+numel(consel), :) = concols(consel, :);

% update UI by re-setting GLM
ch.CMFig.Pointer = cfp;
ne_cm_setglm;
