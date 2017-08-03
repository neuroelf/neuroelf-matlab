% FUNCTION ne_vmp_createsmp: create SMP from current VMP
function varargout = ne_vmp_createsmp(varargin)

% Version:  v1.1
% Build:    16052718
% Date:     May-27 2016, 6:36 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2014, 2015, 2016, Jochen Weber
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
ch = ne_gcfg.h;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% only valid if single VMP and at least one SRF in scenery is selected
if nargin < 3 || numel(varargin{3}) ~= 1 || ~isxff(varargin{3}, {'cmp', 'hdr', 'head', 'vmp'})
    stvar = ne_gcfg.fcfg.StatsVar;
    if numel(stvar) ~= 1 || ~isxff(stvar, true)
        return;
    end
else
    stvar = varargin{3};
end
if nargin < 4 || ~ischar(varargin{4}) || ~strcmpi(varargin{4}(:)', 'cursel')
    mapsel = 1:numel(stvar.MapNames);
else
    mapsel = ne_gcfg.fcfg.StatsVarIdx(:)';
end
scene = ne_gcfg.h.Scenery.UserData;
sceni = ne_gcfg.h.Scenery.Value;
if isempty(scene) || isempty(sceni)

    % get filenames of available files
    fnames = scene(:, end);
    for fc = 1:numel(fnames)
        fnames{fc} = fnames{fc}.FilenameOnDisk;
    end

    % no SRF loaded, SRF files exist?
    csrfs = findfiles(neuroelf_path('colin'), 'colin_midres_*_SPH40k_ICBMnorm.srf', 'depth=1');
    csrfs(~cellfun('isempty', regexp(csrfs, '_BH_'))) = [];
    if numel(csrfs) ~= 2
        csrfs = findfiles(neuroelf_path('colin'), 'colin_midres_*_SPH40k.srf', 'depth=1');
        csrfs(~cellfun('isempty', regexp(csrfs, '_BH_'))) = [];
    end
    if numel(csrfs) ~= 2
        uiwait(warndlg('No surfaces loaded/selected and Colin 40k meshes not available.', ...
            'NeuroElf - info', 'modal'));
        return;
    end

    % load Colin SRFs
    loadcs = questdlg('No surfaces loaded/selected. Load LH/RH Colin 40k meshes?', ...
        'NeuroElf - request', 'Yes', 'No', 'Yes');
    if ~ischar(loadcs) || ~strcmpi(loadcs(:)', 'yes')
        return;
    end

    % load surfaces
    try
        cpage = ne_gcfg.fcfg.page;
        lsrfi = findfirst(strcmpi(fnames, csrfs{1}));
        if isempty(lsrfi)
            lsrf = xff(csrfs{1});
            ne_openfile(0, 0, lsrf, true);
            lsrfi = size(ne_gcfg.h.Scenery.UserData, 1);
        end
        rsrfi = findfirst(strcmpi(fnames, csrfs{2}));
        if isempty(rsrfi)
            rsrf = xff(csrfs{2});
            ne_openfile(0, 0, rsrf, true);
            rsrfi = size(ne_gcfg.h.Scenery.UserData, 1);
        end
        ne_showpage(0, 0, cpage);
        scene = ne_gcfg.h.Scenery.UserData;
        sceni = [lsrfi; rsrfi];
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
        uiwait(errordlg('Error loading/browsing Colin surfaces.', ...
            'NeuroElf - error', 'modal'));
        return;
    end
end

% set pointer and don't allow further callbacks
ne_gcfg.c.incb = true;
ch.MainFig.Pointer = 'watch';
drawnow;

% get VMP filename
[stpath, stname] = fileparts(stvar.FilenameOnDisk);
if isempty(stname)
    stname = sprintf('unsaved.%s', stvar.Filetype);
end

% iterate over SRFs (active in scenery)
try
    cini = ne_gcfg.c.ini;
    copt = struct( ...
        'interp', cini.VMP.CreateSMPInterp, ...
        'ipfrom', cini.VMP.CreateSMPIpFrom, ...
        'ipstep', cini.VMP.CreateSMPIpStep, ...
        'ipto',   cini.VMP.CreateSMPIpTo, ...
        'mapsel', -1, ...
        'method', cini.VMP.CreateSMPMethod);
    for sc = 1:numel(sceni)

        % activate surface
        ne_setcsrf(0, 0, sceni(sc));

        % get filename
        srf = scene{sceni(sc), 4};
        [srpath, srname] = fileparts(srf.FilenameOnDisk);
        if isempty(srname)
            srname = makelabel(scene{sceni(sc), 1});
        end

        % echo
        copt.mapsel = mapsel;
        if ne_gcfg.c.echo
            ne_echo('vmp', 'CreateSMP', scene{sceni(sc), 4}, copt);
        end

        % sample SMP(s)
        if ~isxff(stvar, 'vtc')
            smp = stvar.CreateSMP(scene{sceni(sc), 4}, copt);
            smpext = '.smp';
        else
            smp = stvar.CreateMTC(scene{sceni(sc), 4}, copt);
            smpext = '.mtc';
        end

        % attempt to save
        try
            smp.SaveAs([stname '_SW_' srname smpext]);
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end

        % add to interface
        ne_gcfg.c.incb = false;
        ne_openfile(0, 0, smp, true);
        if nargout >= sc
            varargout{sc} = smp;
        end

        % update stats
        ne_setcsrf(0, 0, sceni(sc));
        ne_gcfg.c.incb = true;
    end
catch ne_eo;
    uiwait(warndlg(['An error occurred: ' ne_eo.message], 'NeuroElf GUI - error', 'modal'));
    return;
end

% return things to normal
ch.MainFig.Pointer = 'arrow';
ne_gcfg.c.incb = false;

% show third page
ne_showpage(0, 0, 3);
