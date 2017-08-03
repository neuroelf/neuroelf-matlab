% FUNCTION ne_cm_addcon: add a contrast
function ne_cm_addcon(varargin)

% Version:  v1.1
% Build:    16040515
% Date:     Apr-05 2016, 3:32 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2016, Jochen Weber
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

% global variable and shortcut structs
global ne_gcfg;
cc = ne_gcfg.fcfg.CM;
ch = ne_gcfg.h.CM.h;

% request new contrast name
cv = [];
if nargin > 2 && ~isempty(varargin{3}) && (ischar(varargin{3}) || ...
   (iscell(varargin{3}) && numel(varargin{3}) == 1 && ischar(varargin{3}{1}) && ~isempty(varargin{3}{1})))
    newcon = varargin{3};
    if ischar(newcon)
        newcon = {newcon(:)'};
    else
        newcon{1} = newcon{1}(:)';
    end
else
    newcon = inputdlg({['Please enter the contrast''s name:' char(10) char(10) ...
        'To use regular expressions against the predictor names, ' ...
        'begin with a $ character, like in $Negative > Neutral']}, ...
        'NeuroElf GUI - input', 1, {'contrast'});
end
if isequal(newcon, 0) || ...
    isempty(newcon)
    return;
end
if iscell(newcon)
    newcon = newcon{1};
end
newcon = ddeblank(newcon);
if isempty(newcon)
    return;
end

% parse contrast?
if newcon(1) == '$'

    % remove $
    newcon = newcon(2:end);
    if isempty(newcon)
        return;
    end

    % > or -
    gt = find(newcon == '>' | newcon == '-');
    if numel(gt) == 1 && gt > 1 && gt < numel(newcon)
        negpart = ddeblank(newcon(gt+1:end));
        newcon = ddeblank(newcon(1:gt-1));
        conname = sprintf('%s > %s', newcon, negpart);
    else
        negpart = '';
        conname = newcon;
    end

    % name
    conname(conname == '$' | conname == '^') = [];
    conname = strrep(strrep(strrep(conname, '.+', '+'), '.*', '*'), '?', '');

    % parse
    prednames = ch.Predictors.String;
    if ~iscell(prednames)
        prednames = cellstr(prednames);
    end
    if ~isempty(strfind(newcon, ' +'))
        newcon = ddeblank(splittocell(newcon, ' +'));
    else
        newcon = {newcon};
    end
    ppredi = false(size(prednames));
    for pc = 1:numel(newcon)
        ppredi = ppredi | (~cellfun('isempty', regexpi(prednames, newcon{pc})));
    end
    if ~isempty(negpart)
        if ~isempty(strfind(negpart, ' +'))
            newcon = ddeblank(splittocell(negpart, ' +'));
        else
            newcon = {negpart};
        end
    else
        newcon = {};
    end
    npredi = false(size(prednames));
    for pc = 1:numel(newcon)
        npredi = npredi | (~cellfun('isempty', regexpi(prednames, newcon{pc})));
    end
    
    % no match?
    if ~any(ppredi) && ~any(npredi)
        return;
    end

    % compute weights
    cv = zeros(size(prednames));
    if any(ppredi)
        cv(ppredi) = 1 / sum(ppredi);
    end
    if any(npredi)
        cv(npredi) = -1 / sum(npredi);
    end

    % replace with better name
    if nargin < 4 || ~ischar(varargin{4}) || isempty(varargin{4})
        newcon = conname;
    else
        newcon = varargin{4}(:)';
    end
end

% put into list of contrasts
cc.cons{end + 1, 1} = newcon;

% already contrasts configured
if ~isempty(cv)
    
    % set new weights to 0
    cc.cons{end, 2} = cv(:);

% fresh contrast
elseif size(cc.cons, 1) > 1

    % set new weights to 0
    cc.cons{end, 2} = zeros(numel(cc.preds), 1);

% no contrasts yet -> take GUI weights
else

    % get current weights
    cc.cons{end, 2} = ne_cm_getweights;
end

% set weights
ne_cm_setweights(cc.cons{end, 2});

% then update dropdown
ch.Contrasts.String = cc.cons(:, 1);
ch.Contrasts.Value = size(cc.cons, 1);
ne_gcfg.h.CM.CMFig.SetGroupEnabled('HasCons', 'on');

% and current GLM
ne_gcfg.fcfg.CM = cc;
ne_cm_updatertv;
ne_cm_updateuis(0, 0, cc.glm);
