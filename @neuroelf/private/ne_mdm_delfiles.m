% FUNCTION ne_mdm_delfiles: delete selected entries from lists
function ne_mdm_delfiles(varargin)

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

% get selection of lists
ls = ch.FuncFiles.Value;
if isempty(ls)
    return;
end

% remove selection from lists
if numel(ch.MParFiles.String) == numel(ch.FuncFiles.String) && ...
   (numel(ch.MParFiles.String) ~= 1 || ...
    ~strcmpi(ch.MParFiles.String{1}, '<no files selected>'))
    ch.MParFiles.String(ls) = [];
end
hFig.SetGroupEnabled('SSelect', 'off');
hFig.SetGroupEnabled('MSelect', 'off');
ch.MParFiles.Value = [];
ch.MParFiles.ListboxTop = min(numel(ch.MParFiles.String), ch.MParFiles.ListboxTop);
ch.FuncFiles.String(ls) = [];
ch.FuncFiles.Value = [];
ch.FuncFiles.ListboxTop = min(numel(ch.FuncFiles.String), ch.FuncFiles.ListboxTop);
ch.DsgnFiles.String(ls) = [];
ch.DsgnFiles.Value = [];
ch.DsgnFiles.ListboxTop = ch.FuncFiles.ListboxTop;

% all empty
if isempty(ch.FuncFiles.String)

    % allow motion parameters to be added!
    ch.UseMotParms.Enable = 'on';
end
