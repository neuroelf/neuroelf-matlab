% FUNCTION ne_mdmvoicondavgrsz: resize MDM VOI condition average window
function ne_mdmvoicondavgrsz(varargin)

% Version:  v0.9c
% Build:    11042919
% Date:     Apr-29 2011, 8:11 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, Jochen Weber
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

% check argument
if nargin < 3 || ...
   ~ischar(varargin{3}) || ...
   ~isfield(ne_gcfg.cc, varargin{3}(:)')
    return;
end

% get figure
cc = ne_gcfg.cc.(varargin{3}(:)');
hFig = cc.Satellite;

% get current position
fsize = hFig.Position(3:4);

% flip OVis flag?
if nargin > 3 && ...
    ischar(varargin{4}) && ...
    strcmpi(varargin{4}(:)', 'optvis')
    ne_gcfg.cc.(varargin{3}(:)').Config.optvis = ...
        (ne_gcfg.cc.(varargin{3}(:)').Tags.OptVis.Value > 0);
    cc = ne_gcfg.cc.(varargin{3}(:)');
    if cc.Config.optvis
        hFig.SetGroupVisible('Opt', 'on');
        cc.Tags.OptVis.Position(1:2) = cc.Config.optvispos;
    else
        hFig.SetGroupVisible('Opt', 'off');
        cc.Tags.OptVis.Position(1:2) = [2, 2];
    end
end

% compute and set new axes position
if cc.Config.optvis
    apos = cc.Config.axpos;
    apos(3:4) = max(1, fsize - (apos(1:2) + [8, 26]));
else
    apos = max(1, [36, 26, (fsize - [44, 46])]);
end
cc.Tags.Axes.Position = apos;
