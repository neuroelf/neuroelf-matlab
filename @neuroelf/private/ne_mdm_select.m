% FUNCTION ne_mdm_select: select from the files
function ne_mdm_select(varargin)

% Version:  v0.9b
% Build:    11111510
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
ch = ne_gcfg.h.MDM.h;

% check arguments
if nargin < 3 || ...
   ~ischar(varargin{3}) || ...
    isempty(varargin{3}) || ...
   ~any(strcmpi(varargin{3}, {'func', 'design', 'mparam'}))
    return;
end

% userps?
userps = strcmpi(ch.MParFiles.Enable, 'on');

% reflect selection
switch lower(varargin{3}(1))
    case {'d'}
        fsel = ch.DsgnFiles.Value;
        ch.FuncFiles.Value = fsel;
        if ~isempty(fsel)
            ch.DsgnFiles.ListboxTop = min(fsel(1), ch.DsgnFiles.ListboxTop);
        end
        ch.FuncFiles.ListboxTop = ch.DsgnFiles.ListboxTop;
        if userps
            ch.MParFiles.Value = fsel;
            ch.MParFiles.ListboxTop = ch.FuncFiles.ListboxTop;
        end
    case {'f'}
        fsel = ch.FuncFiles.Value;
        ch.DsgnFiles.Value = fsel;
        if ~isempty(fsel)
            ch.FuncFiles.ListboxTop = min(fsel(1), ch.FuncFiles.ListboxTop);
        end
        ch.DsgnFiles.ListboxTop = ch.FuncFiles.ListboxTop;
        if userps
            ch.MParFiles.Value = fsel;
            ch.MParFiles.ListboxTop = ch.FuncFiles.ListboxTop;
        end
    case {'m'}
        fsel = ch.MParFiles.Value;
        ch.DsgnFiles.Value = fsel;
        ch.FuncFiles.Value = fsel;
        if ~isempty(fsel)
            ch.MParFiles.ListboxTop = min(fsel(1), ch.MParFiles.ListboxTop);
        end
        ch.FuncFiles.ListboxTop = ch.MParFiles.ListboxTop;
        ch.DsgnFiles.ListboxTop = ch.MParFiles.ListboxTop;
end

% allow QA/QC button
if numel(fsel) == 1
    hFig.SetGroupEnabled('SSelect', 'on');
else
    hFig.SetGroupEnabled('SSelect', 'off');
end

% allow MPPlot button
if userps && ...
    numel(fsel) == 1
    hFig.SetGroupEnabled('MSelect', 'on');
else
    hFig.SetGroupEnabled('MSelect', 'off');
end
