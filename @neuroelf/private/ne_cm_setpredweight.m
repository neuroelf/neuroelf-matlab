% FUNCTION ne_cm_setpredweight: set a given weight to selected predictors
function ne_cm_setpredweight(varargin)

% Version:  v1.1
% Build:    16032913
% Date:     Mar-29 2016, 1:40 PM EST
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
ch = ne_gcfg.h.CM.h;
sel = ch.PredWeights.Value;
if nargin < 3 || (isempty(sel) && ~strcmpi(varargin{3}(:)', 'select'))
    return;
end
if strcmpi(varargin{3}(:)', 'select')
    ccon = ch.Contrasts.String;
    if ~iscell(ccon)
        ccon = cellstr(ccon);
    end
    ccon = ccon{ch.Contrasts.Value};
    if isempty(regexpi(ccon, '^[a-z0-9_]+$'))
        ccon = '  ^$';
    end
    csel = inputdlg({'Regular expression match to select conditions:'}, ...
        'NeuroElf - input', 1, {ccon});
    if ~iscell(csel) || numel(csel) ~= 1 || ~ischar(csel{1}) || isempty(ddeblank(csel{1}))
        return;
    end
    cnames = ch.Predictors.String;
    if ~iscell(cnames)
        cnames = cellstr(cnames);
    end
    csel = find(~cellfun('isempty', regexpi(cnames, ddeblank(csel{1}(:)'))));
    if isempty(csel)
        return;
    end
    cdsel = setdiff(csel(:), sel(:));
    if isempty(cdsel)
        newsel = setdiff(sel(:), csel(:));
    else
        newsel = union(sel(:), csel(:));        
    end
    ch.Predictors.Value = newsel;
    ch.PredWeights.Value = newsel;
    return;
end
w = ch.PredWeights.String;
if ~iscell(w)
    w = cellstr(w);
end
for sc = sel(:)'
    if ~isempty(varargin{3})
        w{sc} = varargin{3};
    else
        selval = inputdlg({'Specify weight value for selected conditions:'}, ...
            'NeuroElf - user input', 1, {'  1'});
        if isempty(selval)
            return;
        end
        try
            selval = str2double(selval);
            if numel(selval) ~= 1|| isinf(selval) || isnan(selval)
                try
                    selval = evalin('base', ['[' selval ']']);
                    if isa(selval, 'double') && numel(selval) == numel(sel) && ...
                       ~any(isinf(selval(:)) | isnan(selval(:))) && sc == sel(1)
                        for scc = sel(:)'
                            w{scc} = sprintf('%g', selval(scc));
                        end
                        ch.PredWeights.String = w;
                        if ~isempty(cc.cons)
                            ne_gcfg.fcfg.CM.cons{ch.Contrasts.Value, 2} = ne_cm_getweights;
                            ne_gcfg.fcfg.CM.glm.RunTimeVars.Contrasts = ne_gcfg.fcfg.CM.cons;
                            ne_cm_updateuis(0, 0, cc.glm);
                        end
                        return;
                    end
                catch ne_eo;
                    rethrow(ne_eo);
                end
                error('BAD_VALUE');
            end
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            return;
        end
        w{sc} = sprintf('%g', selval);
    end
end
ch.PredWeights.String = w;

% update config?
if ~isempty(cc.cons)

    % get weights and set in config
    ne_gcfg.fcfg.CM.cons{ch.Contrasts.Value, 2} = ne_cm_getweights;
    ne_gcfg.fcfg.CM.glm.RunTimeVars.Contrasts = ne_gcfg.fcfg.CM.cons;
    ne_cm_updateuis(0, 0, cc.glm);
end
