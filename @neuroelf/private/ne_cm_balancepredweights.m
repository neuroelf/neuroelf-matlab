% FUNCTION ne_cm_balancepredweights: try to auto-balance contrast so sum:=0
function ne_cm_balancepredweights(varargin)

% Version:  v1.1
% Build:    16022911
% Date:     Feb-29 2016, 11:42 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2014, 2016, Jochen Weber
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
w = ch.PredWeights.String;
if ~iscell(w)
    w = cellstr(w);
end
nw = eval(['[' gluetostringc(w, ';') ']']);
if all(nw == 0)
    return;
end
ppos = (nw > 0);
npos = (nw < 0);
spos = sum(nw(ppos));
nw(ppos) = (1 / spos) .* nw(ppos);
if any(npos)
    sneg = -sum(nw(npos));
    nw(npos) = (1 / sneg) .* nw(npos);
end
bpos = ppos | npos;
for wc = 1:numel(w);
    if bpos(wc)
        w{wc} = sprintf('%g', nw(wc));
    else
        w{wc} = '0';
    end
end
ch.PredWeights.String = w(:);
nw = eval(['[' gluetostringc(w, ';') ']']);

% make sure sums are up to 1, if rounding error occurred
spos = sum(nw(ppos));
if spos ~= 1
    ppos = find(ppos);
    w{ppos(1)} = sprintf('%g', nw(ppos(1)) + 1 - spos);
    ch.PredWeights.String = w(:);
end
if any(npos)
    sneg = -sum(nw(npos));
    if sneg ~= 1
        npos = find(npos);
        w{npos(1)} = sprintf('%g', nw(npos(1)) - (1 - spos));
        ch.PredWeights.String = w(:);
    end
end

% update config?
if ~isempty(cc.cons)

    % get weights and set in config
    ne_gcfg.fcfg.CM.cons{ch.Contrasts.Value, 2} = ne_cm_getweights;
    ne_cm_updatertv;
    ne_cm_updateuis(0, 0, cc.glm);
end
