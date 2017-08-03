function ne_maskstatswithvmr(varargin)
% ne_maskstatswithvmr  - mask the current statsvar with current VMR
%
% FORMAT:       ne_maskstatswithvmr(SRC, EVT [, vtcmode])
%
% Input fields:
%
%       SRC, EVT    Matlab handle callback inputs (discarded)
%       vtcmode     if set to 'vtc' use VTC in ne_gcfg.fcfg.SliceVar
%
% Example:
%
%   neuroelf_gui('openfile', 'results.vmp');
%   vmpcopy = neuroelf_gui('clonefile', 'StatsVar');
%   neuroelf_gui('maskstatswithvmr');

% Version:  v1.0
% Build:    14091315
% Date:     Sep-13 2014, 3:08 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2014, Jochen Weber
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
cc = ne_gcfg.fcfg;

% only valid for certain combinations
if nargin < 3 || ...
   ~ischar(varargin{3}) || ...
   ~strcmpi(varargin{3}(:)', 'vtc')
    if numel(cc.SliceVar) ~= 1 || ...
       ~isxff(cc.SliceVar, 'vmr') || ...
        numel(cc.StatsVar) ~= 1 || ...
       ~isxff(cc.StatsVar, {'glm', 'vmp'}) || ...
       (isxff(cc.StatsVar, 'glm') && ...
        cc.StatsVar.ProjectType ~= 1)
        return;
    end
    stats = true;
else
    if numel(cc.SliceVar) ~= 1 || ...
       ~isxff(cc.SliceVar, 'vtc')
        return;
    end
    stats = false;

    % get list of VMRs
    xr = xff;
    vmrs = xr.Documents('vmr');
    if isempty(vmrs)
        return;
    elseif numel(vmrs) > 1
        vmrf = vmrs;
        for vc = 1:numel(vmrf)
            vmrf{vc} = xr.Document(vmrs{vc}).FilenameOnDisk;
            if isempty(vmrf{vc})
                vmrf{vc} = sprintf('<unsaved.vmr> %d', vc);
            end
        end
        [vmrsel, vmrsok] = listdlg( ...
            'PromptString',  'Please select the VMR to mask with...', ...
            'ListString',    vmrf(:), ...
            'SelectionMode', 'single', ...
            'InitialValue',  1, ...
            'ListSize',      [480, 120]);
        if ~isequal(vmrsok, 1) || ...
            isempty(vmrsel)
            return;
        end
        vmr = xr.Document(vmrs{vmrsel});
    else
        vmr = xr.Document(vmrs{1});
    end
end

% get variable to mask with
if stats
    svar = cc.StatsVar;
    mvar = cc.SliceVar;
else
    svar = cc.SliceVar;
    mvar = vmr;
end

% see if data is transio
svartype = lower(svar.Filetype);
tioobj = false;
switch (svartype)
    case {'glm'}
        if svar.ProjectTypeRFX > 0
            if istransio(svar.GLMData.Subject(1).BetaMaps)
                tioobj = true;
            end
        else
            if istransio(svar.GLMData.BetaMaps)
                tioobj = true;
            end
        end
    case {'vmp'}
        if istransio(svar.Map(1).VMPData)
            tioobj = true;
        end
    case {'vtc'}
        if istransio(svar.VTCData)
            tioobj = true;
        end
end

% request to either load data or save under different name first
if tioobj
    action = questdlg( ...
        ['The stats object''s data is only on disk (and access via transio). ' ...
         'Would you like to load the data or save under a different name first?'], ...
        'NeuroElf - user request', 'Load into memory', 'SaveAs first', 'Cancel', 'Load into memory');
    if isempty(action) || ...
       ~ischar(action) || ...
        any('cC' == action(1))
        return;
    end
    action = lower(action(1));
    if action == 's'
        svarfile = svar.FilenameOnDisk;
        try
            svar.SaveAs;
        catch ne_eo;
            uiwait(warndlg(['Saving failed: ' ne_eo.message], ...
                'NeuroElf - error', 'modal'));
            return;
        end
        if strcmpi(svarfile, svar.FilenameOnDisk)
            uiwait(warndlg('Saved under the same name, not masking.', ...
                'NeuroElf - warning', 'modal'));
            return;
        end
        if svartype(1) == 'g'
            action = questdlg('Preload data (faster masking performance)?', ...
                'NeuroElf - user request', 'Yes', 'No', 'Yes');
            if ischar(action) && ...
               ~isempty(action) && ...
                any('yY' == action(1))
                try
                    svar.LoadTransIOData;
                catch ne_eo;
                    uiwait(warndlg(['Error loading transio data: ' ne_eo.message], ...
                        'NeuroElf - error', 'modal'));
                    return;
                end
            end
        end
    else
        try
            svar.LoadTransIOData;
        catch ne_eo;
            uiwait(warndlg(['Error loading transio data: ' ne_eo.message], ...
                'NeuroElf - error', 'modal'));
            return;
        end
    end
end

% block further calls
ne_gcfg.c.incb = true;
ne_gcfg.h.MainFig.Pointer = 'watch';

% try call
try
    svar.MaskWithVMR(mvar);
    if ne_gcfg.c.echo
        ne_echo(svar.Filetype, 'MaskWithVMR', mvar);
    end
    if isxff(cc.SliceVar, true) && ...
        isfield(handles(cc.SliceVar), 'RenderSVol')
        cc.SliceVar.DeleteHandle('RenderSVol');
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
end

% re-allow calls
ne_gcfg.c.incb = false;
ne_gcfg.h.MainFig.Pointer = 'arrow';

% update display
if cc.page < 3
    ne_setslicepos;
elseif cc.page > 3
    ne_render_setview;
end
