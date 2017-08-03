% FUNCTION ne_cm_loadcons: load contrasts
function ne_cm_loadcons(varargin)

% Version:  v1.1
% Build:    16022911
% Date:     Feb-29 2016, 11:36 AM EST
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

% ask for filename
[ctrfile, ctrpath] = uigetfile( ...
    {'*.ctr',       'Contrasts files (*.ctr)'}, ...
     'Please select the file containing the contrasts...', ...
     'MultiSelect', 'off');
if isequal(ctrfile, 0) || ...
    isequal(ctrpath, 0) || ...
    isempty(ctrfile)
    return;
end
ctrfile = [ctrpath ctrfile];

% try reading contrasts
try
    ctr = xff(ctrfile);
    conn = ctr.ContrastNames(:);
    cons = ctr.ContrastValues';
    ctr.ClearObject;
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    uiwait(warndlg('Invalid contrasts file.', 'NeuroElf GUI - Info', 'modal'));
    return;
end

% get config
cc = ne_gcfg.fcfg.CM;
ch = ne_gcfg.h.CM.h;

% empty contrast file
if isempty(conn) || ...
    isempty(cons)
    uiwait(warndlg('Empty contrasts file.', 'NeuroElf GUI - Info', 'modal'));
    return;
end
ne_gcfg.h.CM.CMFig.SetGroupEnabled('HasCons', 'on');

% set names
cc.cons = cell(numel(conn), 2);
cc.cons(:, 1) = conn(:);

% non RFX
if cc.glm.ProjectTypeRFX == 0

    % number must match one-on-one
    if size(cons, 1) ~= numel(cc.preds)
        uiwait(warndlg('Number of weights mismatch.', 'NeuroElf GUI - Info', 'modal'));
        return;
    end

    % so, store values
    for c = 1:numel(conn)
        cc.cons{c, 2} = cons(:, c);
    end

% for RFX
else

    % numbers must match more complicated
    if size(cons, 1) ~= ((numel(cc.preds) + 1) * cc.nsubs)
        uiwait(warndlg('Number of weights mismatch.', 'NeuroElf GUI - Info', 'modal'));
        return;
    end

    % otherwise equally store (from first subject)
    for c = 1:numel(conn)
        cc.cons{c, 2} = cons(1:numel(cc.preds), c);
    end
end

% set name in control
ch.Contrasts.String = cc.cons(:, 1);
ch.Contrasts.Value = 1;

% and weights
ne_cm_setweights(cc.cons{1, 2});

% and write back in control
ne_gcfg.fcfg.CM = cc;

% and update GLM RunTimeVars
ne_cm_updatertv;
ne_cm_updateuis(0, 0, cc.glm);
