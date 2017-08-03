% FUNCTION ne_cm_addcovs: add covariates
function ne_cm_addcovs(varargin)

% Version:  v0.9d
% Build:    14052913
% Date:     May-29 2014, 1:06 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2014, Jochen Weber
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
ch = ne_gcfg.h.CM.h;

% request name for covariate
newcov = inputdlg({'Please enter the covariate''s name:'}, ...
    'NeuroElf GUI - input', 1, {'covariate'});
if isequal(newcov, 0) || ...
    isempty(newcov)
    return;
end
if iscell(newcov)
    newcov = newcov{1};
end

% request value
newcovval = inputdlg({'Please enter covariate value/variable:'}, ...
    'NeuroElf GUI - input', 1, {''});
if isequal(newcovval, 0) || ...
    isempty(newcovval)
    return;
end
if iscell(newcovval)
    newcovval = newcovval{1};
end

% use evalin('base') to get to contents
try
    newcovval = evalin('base', newcovval);
    if ~isnumeric(newcovval)
        error('MUST_BE_NUMERIC');
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    uiwait(warndlg('Invalid covariate contents/variable.', ...
        'NeuroElf GUI - Info', 'modal'));
    return;
end

% dimensions
if ndims(newcovval) ~= 2 || ...
   ~any(size(newcovval) == cc.nsubs)
    uiwait(warndlg('Covariate must match number of subjects.', ...
        'NeuroElf GUI - Info', 'modal'));
    return;
end

% transpose?
if size(newcovval, 1) ~= cc.nsubs
    newcovval = newcovval';
end

% create name
cvname = cell(1, size(newcovval, 2));
if numel(cvname) > 1
    for cvc = 1:numel(cvname)
        cvname{cvc} = sprintf('%s - %d', newcov, cvc);
    end
else
    cvname{1} = newcov;
end

% add to config
cvc = size(cc.covs, 1);
cc.covs(cvc+1:cvc+numel(cvname), 1) = cvname(:);
for ncvc = 1:numel(cvname)
    cc.covs{cvc + ncvc, 2} = newcovval(:, ncvc);
end

% then set in handle
ch.Covs.String = cc.covs(:, 1);

% and update global array
ne_gcfg.fcfg.CM = cc;

% update in GLM
ne_cm_updatertv;
