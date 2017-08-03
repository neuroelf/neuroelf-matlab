% FUNCTION ne_mdm_closeui: close mdm config dialog
function ne_mdm_closeui(varargin)

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
hFig = ne_gcfg.h.MDM.MDMFig;
cc = ne_gcfg.fcfg.MDM;

% unload MDM
if isfield(cc, 'mdm') && ...
    numel(cc.mdm) == 1 && ...
    isxff(cc.mdm, 'mdm')
    cc.mdm.ClearObject;
end

% clear config
ne_gcfg.fcfg.MDM = [];

% update last known position
ne_gcfg.c.ini.Children.MDMPosition = hFig.Position(1:2);

% delete figure and remove from global struct
hFig.Delete;
ne_gcfg.h.MDM = [];

% release contrast manager call
ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'mdm_open')) = [];
