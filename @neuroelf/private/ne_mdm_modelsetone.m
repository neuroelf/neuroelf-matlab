% FUNCTION ne_mdm_modelsetone: set one modeling way
function ne_mdm_modelsetone(varargin)

% Version:  v0.9c
% Build:    12042611
% Date:     Nov-15 2011, 10:59 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, Jochen Weber
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
if nargin < 3 || ...
   ~isa(varargin{3}, 'double') || ...
    numel(varargin{3}) ~= 1
    return;
end
cf = ne_gcfg.h.MDM.MDMFig;
ch = ne_gcfg.h.MDM.h;

% set one button in radio group
cf.RadioGroupSetOne('Model', varargin{3});

% for random effects and some transformation
if ch.ModelRFX.Value > 0

    % enable VWeight-ing
    ch.VWeight.Enable = 'on';
    ch.VWeight.Value = 1;

    % also enable "from FFX"?
    if ch.Trans0.Value == 0
        ch.CombineFFX.Enable = 'on';
    end

else

    % disable VWeight-ing / "from FFX"
    ch.VWeight.Enable = 'off';
    ch.VWeight.Value = 0;
    ch.CombineFFX.Enable = 'off';
    ch.CombineFFX.Value = 0;
end
