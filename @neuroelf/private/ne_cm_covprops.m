% FUNCTION ne_cm_covprops: update covariate properties
function varargout = ne_cm_covprops(varargin)

% Version:  v1.0
% Build:    16011015
% Date:     Jan-10 2016, 3:23 PM EST
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
cc = ne_gcfg.fcfg.CM;
glm = cc.glm;
ch = ne_gcfg.h.CM.h;
hFig = ne_gcfg.h.CM.CMFig;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% do nothing if multiple covariates selected
covi = ch.Covs.Value;
if numel(covi) ~= 1
    return;
end

% request updated settings
covname = ch.Covs.String{covi};
covset = {covname, '  n', '  n', '', 'new covariate'};
covset = inputdlg({'Covariate name:', 'Set 0 values to NaN? (y/n)', ...
    'Remove mean? (y/n)', 'Compute new covariate from expression:', ...
    'And store as:'}, ...
    'NeuroElf GUI - Map settings', 1, covset);
if ~iscell(covset) || ...
    numel(covset) ~= 5
    return;
end

% update name
if ~strcmp(ddeblank(covset{1}), covname) && ...
   ~isempty(ddeblank(covset{1}))
    glm.RunTimeVars.CovariatesNames{covi} = ddeblank(covset{1});
    ne_gcfg.fcfg.CM.covs{covi, 1} = covset{1};
    ch.Covs.String{covi} = covset{1};
end

% set 0 values to NaN
cvals = glm.RunTimeVars.CovariatesData(:, covi);
upvals = false;
if ~isempty(covset{2}) && ...
    any(lower(covset{2}) == 'y')
    cvals(cvals == 0) = NaN;
    upvals = true;
end

% remove mean
if ~isempty(covset{3}) && ...
    any(lower(covset{3}) == 'y')
    cvals = cvals - meannoinfnan(cvals);
    upvals;
end

% update
if upvals
    glm.RunTimeVars.CovariatesData(:, covi) = cvals;
    ne_gcfg.fcfg.CM.covs{covi, 2} = cvals;
    cc = ne_gcfg.fcfg.CM;
end

% and compute formula?
if ~isempty(ddeblank(covset{4})) && ...
   ~isempty(ddeblank(covset{5}))
    
    % assign local variables
    for cvc = size(cc.covs, 1):-1:1
        eval(sprintf('%s=cc.covs{cvc, 2};', makelabel(lower(cc.covs{cvc, 1}))));
    end
    
    % try to evaluate formula
    try
        newcvals = eval(lower(ddeblank(covset{4})));
    catch ne_eo;
        uiwait(warndlg(sprintf('Error evaluating expression\n%s:\n\n%s', ...
            ddeblank(covset{4}), ne_eo.message), 'modal'));
        if upvals
            ne_cm_updateuis(0, 0, cc.glm);
        end
        return;
    end
    
    % add to GLM and interface
    if isequal(size(newcvals), size(cvals)) && ...
       ~all(isnan(newcvals)) && ...
       ~any(isinf(newcvals)) && ...
        var(newcvals) > sqrt(eps)
        glm.RunTimeVars.CovariatesNames{end+1} = ddeblank(covset{5});
        glm.RunTimeVars.CovariatesData(:, end+1) = newcvals;
        ch.Covs.String{numel(glm.RunTimeVars.CovariatesNames)} = ddeblank(covset{5});
        ch.Covs.Value = numel(glm.RunTimeVars.CovariatesNames);
        ne_gcfg.fcfg.CM.covs(end+1, :) = ...
            {glm.RunTimeVars.CovariatesNames{end+1}, newcvals};
        upvals = true;
    end
end

% update interfaces
if upvals
    ne_cm_updateuis(0, 0, cc.glm);
end
