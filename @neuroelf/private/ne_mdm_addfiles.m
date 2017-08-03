% FUNCTION ne_mdm_addfiles: manually add one set of files
function ne_mdm_addfiles(varargin)

% Version:  v0.9d
% Build:    14061710
% Date:     Jun-17 2014, 10:04 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2014, Jochen Weber
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

% use RPs ?
userps = strcmpi(ch.MParFiles.Enable, 'on');

% request filenames
[vtcfile, vtcfld] = uigetfile( ...
    {'*.vtc', 'BrainVoyager QX VTC files (*.vtc)'}, ...
    'Please select the VTC file to add...', 'MultiSelect', 'off');
if isequal(vtcfile, 0) || ...
    isequal(vtcfld, 0) || ...
   ~ischar(vtcfile) || ...
   ~ischar(vtcfld) || ...
    isempty(vtcfile) || ...
    isempty(vtcfld)
    return;
end
[desfile, desfld] = uigetfile( ...
    {'*.prt', 'BrainVoyager QX Protocol files (*.prt)'; ...
     '*.rtc', 'BrainVoyager QX Reference time course files (*.rtc)'; ...
     '*.sdm', 'BrainVoyager QX Single-study design matrix files (*.sdm)'}, ...
    'Please select the design file to add...', 'MultiSelect', 'off');
if isequal(desfile, 0) || ...
    isequal(desfld, 0) || ...
   ~ischar(desfile) || ...
   ~ischar(desfld) || ...
    isempty(desfile) || ...
    isempty(desfld)
    return;
end
if userps
    [rpsfile, rpsfld] = uigetfile( ...
        {'*.txt', 'SPM realignment parameter files (*.txt)'; ...
         '*.rtc', 'BrainVoyager QX Reference time course files (*.rtc)'; ...
         '*.sdm', 'BrainVoyager QX Single-study design matrix files (*.sdm)'}, ...
        'Please select the design file to add...', 'MultiSelect', 'off');
    if isequal(rpsfile, 0) || ...
        isequal(rpsfld, 0) || ...
       ~ischar(rpsfile) || ...
       ~ischar(rpsfld) || ...
        isempty(rpsfile) || ...
        isempty(rpsfld)
        return;
    end
end

% make sure this correctly overrides
if strcmpi(ch.FuncFiles.Enable, 'off') || ...
   (numel(ch.FuncFiles.String) == 1 && ...
    strcmpi(ch.FuncFiles.String{1}, '<no files selected>'))
    ch.FuncFiles.String = {};
    ch.FuncFiles.Enable = 'on';
end
if strcmpi(ch.DsgnFiles.Enable, 'off') || ...
   (numel(ch.DsgnFiles.String) == 1 && ...
    strcmpi(ch.DsgnFiles.String{1}, '<no files selected>'))
    ch.DsgnFiles.String = {};
    ch.DsgnFiles.Enable = 'on';
end

% add to end of list and set new index
ch.FuncFiles.String = [ch.FuncFiles.String; {[vtcfld vtcfile]}];
ch.DsgnFiles.String = [ch.DsgnFiles.String; {[desfld desfile]}];
ch.FuncFiles.Value = numel(ch.FuncFiles.String);
ch.DsgnFiles.Value = numel(ch.DsgnFiles.String);
hFig.SetGroupEnabled('SSelect', 'on');
if userps
    if strcmpi(ch.MParFiles.Enable, 'off') || ...
       (numel(ch.MParFiles.String) == 1 && ...
        strcmpi(ch.MParFiles.String{1}, '<no files selected>'))
        ch.MParFiles.String = {};
        ch.MParFiles.Enable = 'on';
        hFig.SetGroupEnabled('MSelect', 'on');
    end
    ch.MParFiles.String = [ch.MParFiles.String; {[rpsfld rpsfile]}];
    ch.MParFiles.Value = numel(ch.MParFiles.String);
else
    ch.UseMotParms.Enable = 'off';
    hFig.SetGroupEnabled('MSelect', 'off');
end
if ~any(cellfun('isempty', regexpi(ch.DsgnFiles.String, '\.prt$')))
    cf.SetGroupEnabled('PRTDsgn', 'on');
else
    cf.SetGroupEnabled('PRTDsgn', 'off');
end
