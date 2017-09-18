% FUNCTION ne_importvmp: call importvmpfromspms
function varargout = ne_importvmp(varargin)

% Version:  v1.1
% Build:    17091810
% Date:     Sep-18 2017, 10:27 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2017, Jochen Weber
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

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% block any callbacks (doesn't update the screen, so this is OK)
if ne_gcfg.c.incb
    return;
end
ne_gcfg.c.incb = true;

% import into current VMP check
if numel(ne_gcfg.fcfg.StatsVar) == 1 && ...
    isxff(ne_gcfg.fcfg.StatsVar, 'vmp')
    ia = questdlg('Import into current VMP?', 'NeuroElf - request', ...
        'Yes', 'No', 'Cancel', 'Yes');
    if isempty(ia) || ~ischar(ia) || strcmpi(ia(:)', 'cancel')
        ne_gcfg.c.incb = false;
        return;
    end
    intocvmp = strcmpi(ia(:)', 'yes');
else
    intocvmp = false;
end

% conversion
try

    % import into current VMP
    if intocvmp

        % create cell with empty argument for clearing in catch block
        vmp = cell(1, 1);

        % select files for import
        if mainver > 6
            msargs = {'MultiSelect', 'on'};
        else
            msargs = {};
        end
        [mapf, mapp] = uigetfile( ...
            {'*.hdr;*.nii;*.nii.gz', 'SPM result map files (*.hdr, *.nii, *.nii.gz)'; ...
             '*.mat', 'SPM result file (*.mat)'}, ...
            'Please select the SPM map files or SPM.mat file to import...', '*.hdr', msargs{:});
        if isequal(mapf, 0) || isequal(mapp, 0)
            ne_gcfg.c.incb = false;
            return;
        end
        if ~iscell(mapf)
            maps = {strrep([strrep(mapp, '\', '/') '/' mapf], '//', '/')};
        else
            maps = mapf;
            for mc = 1:numel(maps)
                maps{mc} = strrep([strrep(mapp, '\', '/') '/' maps{mc}], '//', '/');
            end
        end
        mc = 1;
        maps = maps(:);
        while mc <= numel(maps)
            if ~isempty(regexpi(maps{mc}, '\.mat$'))
                try
                    spmc = load(maps{mc});
                    if ~isfield(spmc, 'SPM')
                        error('NOT_AN_SPM_FILE');
                    end
                    spmp = fileparts(maps{mc});
                    coni = {spmc.SPM.xCon.Vspm};
                    for smc = 1:numel(coni)
                        [conip, conif, conie] = fileparts(coni{smc}.fname);
                        coni{smc} = [spmp '/' conif conie];
                    end
                    maps = [maps(1:mc-1); coni(:); maps(mc+2:end)];
                    mc = mc + numel(coni) - 1;
                catch ne_eo;
                    ne_gcfg.c.lasterr = ne_eo;
                    ne_gcfg.c.incb = false;
                    return;
                end
            end
            mc = mc + 1;
        end

        % call import function
        ne_gcfg.fcfg.StatsVar.ImportSPMMaps(maps, struct( ...
            'interp',  ne_gcfg.c.ini.Statistics.ImportInterpMethod, ...
            'maptype', 'a'));

        % after success, set object
        vmp{1} = ne_gcfg.fcfg.StatsVar;

        % don't alter loaded flag
        loaded = {};

    % create new VMP
    else
        vmp = cell(1, 1);
        vmp{1} = importvmpfromspms([], 'a', [], [], ...
            ne_gcfg.c.ini.Statistics.ImportInterpMethod);

        % if not a VMP, return
        if numel(vmp{1}) ~= 1 || ~isxff(vmp{1}, 'vmp')

            % re-allow callbacks
            ne_gcfg.c.incb = false;
            return;
        end

        % request target file(s)
        [vmpfile, vmppath] = uiputfile( ...
           {'*.vmp', 'BrainVoyager Volume MaP files'}, ...
            'Please specify the target VMP filename...', '');

        % not to be saved ?
        if isequal(vmpfile, 0) || ...
            isequal(vmppath, 0)

            % and re-allow callbacks and return
            ne_gcfg.c.incb = false;

            % then simply add to interface
            ne_openfile(0, 0, vmp{1}, true);
            return;
        end

        % save VMPs
        if iscell(vmppath)
            vmppath = vmppath{1};
        end
        if iscell(vmpfile)
            vmpfile = vmpfile{1};
        end
        vmp{1}.SaveAs(strrep( ...
            [strrep(vmppath, '\', '/') '/' vmpfile], '//', '/'));

        % loaded by function
        loaded = {true};
    end
% on any error
catch ne_eo;

    % clear attempted object
    clearxffobjects(vmp);

    % give error dialog
    errordlg(['Error converting spmX_xxxx.img -> VMP: ' ne_eo.message], ...
        'NeuroElf GUI - error', 'modal');

    % re-allow callbacks
    ne_gcfg.c.incb = false;

    % and get out
    return;
end

% simply re-allow callbacks
ne_gcfg.c.incb = false;

% then open
ne_openfile(0, 0, vmp{1}, loaded{:});
