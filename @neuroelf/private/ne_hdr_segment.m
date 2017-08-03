% FUNCTION ne_hdr_segment: apply SPM segmentation to HDR
function ne_hdr_segment(varargin)

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:34 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2014, 2016, Jochen Weber
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

% requires current object to be HDR
if numel(cc.SliceVar) ~= 1 || ...
   ~isxff(cc.SliceVar, 'hdr')
    return;
end

% try using SPM
try
    spmver = 'unknown';
    spmver = lower(spm('ver'));
    if ~any(strcmp(spmver, {'spm5'; 'spm8', 'spm12b', 'spm12'}))
        error( ...
            'neuroelf:UnsupportedSPM', ...
            'Unsupported version of SPM.' ...
        );
    end
catch ne_eo;
    warndlg(sprintf('SPM error (version %s): %s.', spmver, ne_eo.message), ...
        'NeuroElf - warning', 'modal');
    return;
end

% temporary file?
istemp = false;
hdr = cc.SliceVar;
hdr.LoadTransIOData;
hdrfullfile = hdr.FilenameOnDisk;
if isempty(hdrfullfile)
    hdrfullfile = [tempname '.nii'];
    hdr.FileMagic = 'n+1';
    hdr.NIIFileType = 2;
    try
        hdr.SaveAs(hdrfullfile);
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
        warndlg('Error saving HDR object to temp file.', 'NeuroElf - warning', 'modal');
        return;
    end
    istemp = true;
end

% run segmentation
try
    delfiles = lower(questdlg('Delete temporary files (GM/WM/CSF + modulated)?', ...
        'NeuroElf - user request'));
    if delfiles(1) == 'c'
        if istemp
            delete(hdrfullfile);
        end
        return;
    end
    spmx_extract_brain(hdrfullfile, true, delfiles(1) == 'n');
catch ne_eo;
    warndlg(ne_eo.message, 'NeuroElf - warning', 'modal');
    return;
end

% delete temp files we don't want
[hdrpath, hdrfile, hdrext] = fileparts(hdrfullfile);
xhdrfile = [hdrpath '/x' hdrfile hdrext];
wxhdrfile = [hdrpath '/wx' hdrfile hdrext];
segsnfile = [hdrpath '/' hdrfile '_seg_sn.mat'];
seginvsnfile = [hdrpath '/' hdrfile '_seg_inv_sn.mat'];
if istemp
    delete(hdrfullfile);
    delete(wxhdrfile);
    delete(segsnfile);
    delete(seginvsnfile);
end

% load segmentated file
try
    xhdr = xff(xhdrfile);
    xhdr.LoadTransIOData;
catch ne_eo;
    if istemp
        delete(xhdrfile);
    end
    warndlg(ne_eo.message, 'NeuroElf - warning', 'modal');
    return;
end

% temporary file?
if istemp

    % delete file
    delete(xhdrfile);

    % copy orientation (from xhdr to hdr)
    xhdr.Browse;
    hdr.RunTimeVars.Mat44 = xhdr.RunTimeVars.Mat44;
    hdr.RunTimeVars.TrfPlus = xhdr.RunTimeVars.TrfPlus;
    hdr.RunTimeVars.Trf = xhdr.RunTimeVars.Trf;
    hdr.DataHist.NIftI1 = xhdr.DataHist.NIftI1;

    % create copy
    newhdr = hdr.CopyObject;

    % copy data
    newhdr.VoxelData = xhdr.VoxelData;

    % delete old object and rename
    xhdr.UnBrowse;
    xhdr.ClearObject;
    xhdr = newhdr;

% otherwise
else

    % reload from disk (coregistration update)
    hdr.ReloadFromDisk;
    hdr.LoadTransIOData;
end

% get "skull only" data
datares = hdr.CoordinateFrame.Resolution(1:3);
brain = (xhdr.VoxelData > 0);
skull = hdr.VoxelData;
skull(brain) = 0;

% grow brain
head = brain;
for gc = 1:ceil(8 / harmmean(datares))
    head = (smoothdata3(single(head), [2, 2, 2]) > 0.25);
end

% clean skull data
sskull = smoothdata3(single(skull), [2, 2, 2] ./ datares);

% browse both objects (update screen)
hdr.Browse;
xhdr.Browse;
