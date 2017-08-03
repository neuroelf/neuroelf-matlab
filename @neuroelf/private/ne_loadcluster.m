function varargout = ne_loadcluster(varargin)
% ne_loadcluster  - load cluster (VOI) file
%
% FORMAT:       ne_loadcluster(SRC, EVT [, voifile])
%
% Input fields:
%
%       SRC, EVT    Matlab handle callback inputs (discarded)
%       voifile     either filename or VOI object, requested if not given
%
% No output fields.
%
% Example:
%
%     ne_loadcluster(0, 0, [neuroelf_path('files') '/shenparcel/shen_parcels.voi']);

% Version:  v1.1
% Build:    16121617
% Date:     Dec-16 2016, 5:08 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2016, Jochen Weber
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

% no filename/object given
loaded = true;
if nargin < 3 || ((~ischar(varargin{3}) || isempty(varargin{3})) && ...
    (numel(varargin{3}) ~= 1 || ~isxff(varargin{3}, 'voi')))

    % what text
    if ne_gcfg.c.ini.Statistics.ClusterTableAdd
        nftitle = 'Please select a cluster definition file to add...';
    else
        nftitle = 'Please select a cluster definition file to load...';
    end

    % try to get filename
    [newvoifile, nvpath] = uigetfile( ...
        {'*.voi', 'Volume-of-Interest file (*.voi)'; ...
         '*.hdr;*.nii;*.nii.gz', 'Analyze mask/cluster file (*.hdr, *.nii, *.nii.gz)'; ...
         '*.msk', 'BrainVoyager QX mask file (*.msk)'; ...
         '*.txt', 'Cluster definitions (*.txt)'}, nftitle);
    if isequal(nvpath, 0) || isequal(newvoifile, 0) || isempty(newvoifile) || ...
        isempty(regexpi(newvoifile, '\.(hdr|msk|nii|nii\.gz|txt|voi)$'))
        return;
    end
    if isempty(nvpath)
        nvpath = pwd;
    end
    newvoifile = [nvpath filesep newvoifile];
    if ~isempty(regexpi(newvoifile, '\.(hdr|nii)$'))
        ithresh = inputdlg('Threshold value for image:', 'NeuroElf - input', 1, {'1'});
        if ~iscell(ithresh) || numel(ithresh) ~= 1 || ...
            isempty(regexpi(ithresh{1}, '^[0-9\.]+$'))
            return;
        end
        ithresh = str2double(ithresh) - sqrt(eps);
    else
        ithresh = 1 - sqrt(eps);
    end
    try
        if ~strcmpi(newvoifile(end-2:end), 'voi')
            sepclus = questdlg('Import as separate clusters?', ...
                'NeuroElf request', 'Yes', 'No', 'Cancel', 'Yes');
            if ~ischar(sepclus) || ~any(strcmpi(sepclus, {'yes', 'no'}))
                return;
            end
            sepclus = strcmpi(sepclus, 'yes');
            lfile = '';
            lsort = 'size';
            if sepclus
                [lfilename, lfilepath] = uigetfile( ...
                    {'*.txt', 'Label text file (*.txt)'}, ...
                    'Please select a label file... (press cancel for no labels)');
                if isequal(lfilename, 0) || isequal(lfilepath, 0)
                    lfile = '';
                else
                    if isempty(lfilepath)
                        lfilepath = pwd;
                    end
                    lfile = [lfilepath filesep lfilename];
                    lsort = 'value';
                end
            end
            newvoi = xff('new:voi');
            newvoi.ImportClusters(newvoifile, ...
                struct('ithresh', ithresh, 'sepclus', sepclus, 'labels', lfile, 'sort', lsort));
        else
            newvoi = xff(newvoifile);
        end
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
        uiwait(warndlg(['Error loading VOIs from file: ' ne_eo.message], ...
            'NeuroElf GUI - I/O error', 'modal'));
        newvoi = [];
    end
    if numel(newvoi) ~= 1
        return;
    end
    if isxff(newvoi) && ...
       ~isxff(newvoi, 'voi')
        newvoi.ClearObject;
        return;
    end

% filename given
elseif ischar(varargin{3})

    % try loading
    try
        newvoi = cell(1, 1);
        newvoi{1} = xff(varargin{3}(:)');
        if numel(newvoi{1}) ~= 1
            return;
        end
        if isxff(newvoi{1}) && ~isxff(newvoi{1}, 'voi')
            clearxffobjects(newvoi);
            return;
        end
        newvoi = newvoi{1};
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
        clearxffobjects(newvoi);
        return;
    end

% object given
else
    loaded = false;
    newvoi = varargin{3};
end

% echo
if ne_gcfg.c.echo && ~isempty(newvoi.FilenameOnDisk)
    ne_echo({'voi = xff(''%s'');', newvoi.FilenameOnDisk});
end

% content being added to VOIs
if ne_gcfg.c.ini.Statistics.ClusterTableAdd

    % try to add to VOI
    try
        ne_gcfg.voi.VOI = joinstructs(ne_gcfg.voi.VOI(:), newvoi.VOI(:));
        ne_gcfg.voi.NrOfVOIs = numel(ne_gcfg.voi.VOI);
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
        uiwait(warndlg(['Error adding VOIs: ', ne_eo.message], ...
            'Neuroelf - error', 'modal'));
        if loaded
            newvoi.ClearObject;
        end
        return;
    end
    if loaded
        newvoi.ClearObject;
    end

% content replacing VOIs
else

    % clear current VOI
    try
        voih = handles(ne_gcfg.voi);
        if isfield(voih, 'RGBImage') && ...
            numel(voih.RGBImage) == 1 && ...
            isxff(voih.RGBImage, 'hdr')
            voih.RGBImage.ClearObject;
        end
        ne_gcfg.voi.ClearObject;
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end

    % take object if loaded
    if loaded
        ne_gcfg.voi = newvoi;

    % otherwise, make a copy
    else
        ne_gcfg.voi = newvoi.CopyObject;
    end
end

% build new list of names
voi = ne_gcfg.voi;
clnames = voi.VOINames;
for nc = 1:numel(clnames)
    cl = voi.VOI(nc);
    clnames{nc} = sprintf('%s (%d voxels around [%d, %d, %d])', ...
        cl.Name, size(cl.Voxels, 1), round(mean(cl.Voxels, 1)));
end
ch.Clusters.ListboxTop = 1;
ch.Clusters.Value = [];
ch.Clusters.String = clnames;
ch.Clusters.Enable = 'on';
