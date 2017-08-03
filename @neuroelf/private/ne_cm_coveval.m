% FUNCTION ne_cm_coveval: re-code covariates by expression
function varargout = ne_cm_coveval(varargin)

% Version:  v1.1
% Build:    16081714
% Date:     Aug-17 2016, 2:58 PM EST
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

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% no expression given
covdata = glm.RunTimeVars.CovariatesData;
covnames = strrep(ch.Covs.String, ' ', '_');
if isempty(covdata) || isempty(covnames)
    return;
end
if nargin < 4 || ~ischar(varargin{3}) || isempty(varargin{3}) || ...
   ~ischar(varargin{4}) || isempty(varargin{4}) || ~any(varargin{4}(:) == '$')
    xstr = inputdlg({'Name of new (re-coded) covariate:', ...
        'New covariate expression: (use $NAME for existing variables)'}, ...
        'NeuroElf GUI - recode covariate', 1, {'newcov', sprintf('  $%s + 1', covnames{1})});
    if ~iscell(xstr) || numel(xstr) ~= 2 || isempty(ddeblank(xstr{1})) || ...
        isempty(ddeblank(xstr{2})) || ~any(xstr{2} == '$')
        return;
    end
    nstr = ddeblank(xstr{1});
    xstr = ddeblank(xstr{2});
else
    nstr = ddeblank(varargin{3}(:)');
    xstr = ddeblank(varargin{4}(:)');
end

% replace covariate names
while any(xstr == '$')
    xrstr = regexpi(xstr, '\$(\w+)', 'tokens');
    if isempty(xrstr)
        return;
    end
    for rc = 1:numel(xrstr)
        if ~any(strcmpi(covnames, xrstr{rc}{1}))
            uiwait(warndlg(['Covariate ' xrstr{rc}{1} ' not found.'], ...
                'NeuroElf - error', 'modal'));
            return;
        end
        xstr = strrep(xstr, ['$' xrstr{rc}{1}], sprintf('covdata(:,%d)', ...
            findfirst(strcmpi(covnames, xrstr{rc}{1}))));
    end
end

% try to evaluate
try
    newcovdata = eval(['zeros(size(covdata,1),1) + (' xstr ')']);
catch ne_eo;
    uiwait(warndlg(['Error evaluating expression:', char([10, 10]), ...
        ne_eo.message], 'NeuroElf - error', 'modal'));
    return;
end

% update GLM
glm.RunTimeVars.CovariatesNames(end+1) = {nstr};
glm.RunTimeVars.CovariatesData(:, end+1) = newcovdata;
glm.RunTimeVars.CovariatesColors(end+1, :) = ...
    {[32, 32, 255; 0, 160, 128; 128, 192, 0; 224, 144, 0; 255, 64, 0], [NaN, NaN]};

% update interfaces
ne_gcfg.fcfg.CM.covs(end+1, :) = {nstr, newcovdata};
covnames = ch.Covs.String;
covnames{end+1} = nstr;
ch.Covs.String = covnames;
ne_cm_updateuis(0, 0, cc.glm);
