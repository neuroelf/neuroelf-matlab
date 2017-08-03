% FUNCTION ne_cm_savecons: save contrasts
function ne_cm_savecons(varargin)

% Version:  v0.9b
% Build:    10081109
% Date:     Aug-11 2010, 9:00 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, Jochen Weber
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

% compile CTR object
cons = cc.cons;
ctr = xff('new:ctr');
ctr.NrOfContrasts = size(cons, 1);

% for non-RFX
if cc.glm.ProjectTypeRFX == 0

    % set initial values
    ctr.NrOfValues = numel(cc.preds);
    ctr.NrOfContrasts = size(cons, 1);
    ctr.ContrastNames = cons(:, 1)';
    ctr.ContrastValues = zeros(ctr.NrOfContrasts, ctr.NrOfValues);

    % copy over
    for c = 1:ctr.NrOfContrasts
        ctr.ContrastValues(c, :) = cc.cons{c, 2}(:)';
    end

% for RFX
else

    % set initial values
    ctr.NrOfValues = (numel(cc.preds) + 1) * cc.nsubs;
    ctr.NrOfContrasts = size(cons, 1);
    ctr.ContrastNames = cons(:, 1)';
    ctr.ContrastValues = zeros(ctr.NrOfContrasts, ctr.NrOfValues);

    % copy over
    for c = 1:ctr.NrOfContrasts
        ctr.ContrastValues(c, :) = repmat([cc.cons{c, 2}(:)', 0], 1, cc.nsubs);
    end
end

% hand over control to SaveAs
ctr.SaveAs;

% then clear object
ctr.ClearObject;
