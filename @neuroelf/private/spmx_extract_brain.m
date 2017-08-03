function spmx_extract_brain(file, keepnpar, keepfiles)
% spmx_extract_brain  - extract the brain from an anatomical image
%
% FORMAT:       spmx_extract_brain(anatfile [, keepnpar [, keepfiles]])
%
% Input fields:
%
%       anatfile    anatomical image filename
%       keepnpar    keep normalization parameters (default: true)
%       keepfiles   keep additional files (Cx + modulated, false)
%
% No output fields.
%
% Note: the file will have an "x" prepended to the filename.

% Version:  v0.9d
% Build:    14081411
% Date:     Aug-14 2014, 11:00 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2014, Jochen Weber
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

% argument check
if nargin > 0 && ...
    iscell(file) && ...
    numel(file) == 1
    file = file{1};
end
if nargin < 1 || ...
   ~ischar(file) || ...
    isempty(file) || ...
    exist(file, 'file') ~= 2
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument or file not found.' ...
    );
end
if nargin < 2 || ...
   ~islogical(keepnpar) || ...
    numel(keepnpar) ~= 1
    keepnpar = true;
end
if nargin < 3 || ...
   ~islogical(keepfiles) || ...
    numel(keepfiles) ~= 1
    keepfiles = false;
end

% detect SPM version
try
    spmv = str2double(regexprep(lower(spm('ver')), '^spm(\d+).*$', '$1'));
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    error( ...
        'neuroelf:SPMError', ...
        'Missing SPM installation or invalid SPM version detected.' ...
    );
end
if numel(spmv) ~= 1 || ...
    isinf(spmv) || ...
    isnan(spmv) || ...
   ~any(spmv == [5, 8, 12])
    error( ...
        'neuroelf:SPMError', ...
        'Unsupported SPM version detected.' ...
    );
end

% get template filenames
spmp = fileparts(which('spm'));
t1m = [spmp filesep 'templates' filesep 'T1.nii'];
if spmv > 8 && ...
    exist(t1m, 'file') ~= 2
    t1m = [spmp filesep 'canonical' filesep 'avg305T1.nii'];
end
if exist(t1m, 'file') ~= 2
    error( ...
        'neuroelf:SPMError', ...
        'SPM based T1 template file not found.' ...
    );
end
tpm = { ...
    [spmp filesep 'tpm' filesep 'grey.nii']; ...
    [spmp filesep 'tpm' filesep 'white.nii']; ...
    [spmp filesep 'tpm' filesep 'csf.nii']};
if spmv > 8 && ...
    exist(tpm{3}, 'file') ~= 2
    tpm = { ...
        [spmp filesep 'tpm' filesep 'TPM.nii,1']; ...
        [spmp filesep 'tpm' filesep 'TPM.nii,2']; ...
        [spmp filesep 'tpm' filesep 'TPM.nii,3']};
    if exist(tpm{1}(1:end-2), 'file') ~= 2
        error( ...
            'neuroelf:SPMError', ...
            'SPM based TPM files not found.' ...
        );
    end
elseif exist(tpm{1}, 'file') ~= 2 || ...
    exist(tpm{2}, 'file') ~= 2 || ...
    exist(tpm{3}, 'file') ~= 2
    error( ...
        'neuroelf:SPMError', ...
        'SPM based TPM files not found.' ...
    );
end

% if result exists, exit
opath = pwd;
[fpath, fname, fext] = fileparts(file);
if isempty(fpath)
    fpath = opath;
end
tfile = [fpath '/x' fname fext];
if exist(tfile, 'file') > 0
    return;
end

% check extension and replace with hdr for xff access later
xext = fext;
if ~isempty(regexpi(xext, 'img$'))
    xext = regexprep(xext, 'img$', 'hdr', 'preservecase');
end

% load job
if spmv < 8
    error( ...
        'neuroelf:NotYetImplemented', ...
        'Brain extraction only available for SPM8.' ...
    );
    % crjob = load([neuroelf_path('spm') '/spm5_extract.mat']);
else
    crjob = load([neuroelf_path('spm') '/spm8_extract.mat']);
end
crjob = crjob.jobs;

% adapt job
file = [file ',1'];
if spmv < 8
    crjob{1}.spatial{1}.smooth{1}.data{1} = file;
    crjob{2}.spatial{1}.coreg{1}.estimate.ref{1} = [t1m ',1'];
    crjob{3}.spatial{1}.preproc{1}.opts.tpm = tpm;
else
    crjob{1}.spm.spatial.smooth.data{1} = file;
    crjob{2}.spm.spatial.coreg.estimate.other{1} = file;
    crjob{2}.spm.spatial.coreg.estimate.ref{1} = [t1m ',1'];
    crjob{3}.spm.spatial.preproc.data{1} = file;
    crjob{3}.spm.spatial.preproc.opts.tpm = tpm;
    if spmv > 8
        crjob{3}.spm.tools.oldseg = crjob{3}.spm.spatial.preproc;
        crjob{3}.spm = rmfield(crjob{3}.spm, 'spatial');
    end
end

% initialize SPM
if spmv < 8
    spm_defaults;
else
    spm('defaults', 'FMRI');
    spm_jobman('initcfg');
end

% run jobs
cd(fpath);
spm_jobman('run', crjob);

% mid-clean up
delete([fpath '/s' fname fext]);
if ~strcmp(fext, xext)
    try
        delete([fpath '/s' fname xext]);
    catch ne_eo;
        warning( ...
            'neuroelf:DeleteFailed', ...
            'Error deleting IMG file: %s', ...
            ne_eo.message ...
        );
    end
end

% get data
modulated = xff([fpath '/m' fname xext]);
class1 = xff([fpath '/c1' fname xext]);
class2 = xff([fpath '/c2' fname xext]);
class3 = xff([fpath '/c3' fname xext]);

% create mask
class123 = (uint16(class1.VoxelData(:, :, :)) + uint16(class2.VoxelData(:, :, :)) + uint16(class3.VoxelData(:, :, :))) > uint16(127);

% close objects
class1.ClearObject;
class2.ClearObject;
class3.ClearObject;

% smooth and cluster mask
class123 = smoothdata3(double(class123), [2, 2, 2]) > 0.5;
[cs, class123] = clustercoordsc(class123, 2, 8);
csmax = maxpos(cs);
class123 = (class123 == csmax);
[cs, class123] = clustercoordsc(~class123, 2, 8);
csmax = maxpos(cs);
class123 = (class123 ~= csmax);

% load and mask data
modulated.LoadVoxelData;
nmodulated = modulated.CopyObject;
modulated.ClearObject;
nmodulated.VoxelData(~class123) = 0;

% save as brain extracted
nmodulated.SaveAs([fpath '/x' fname xext]);
nmodulated.ClearObject;

% normalize brain
spm_write_sn([fpath '/x' fname fext], [fpath '/' fname '_seg_sn.mat'], ...
    struct('interp', 3, 'wrap', [0, 0, 0], 'vox', [1, 1, 1], ...
    'bb', [-90, -126, -72; 90, 90, 102]));

% clean up...
if ~keepfiles
    delete([fpath '/c1' fname fext]);
    delete([fpath '/c2' fname fext]);
    delete([fpath '/c3' fname fext]);
    delete([fpath '/m' fname fext]);
    if ~strcmp(fext, xext)
        try
            delete([fpath '/c1' fname xext]);
            delete([fpath '/c2' fname xext]);
            delete([fpath '/c3' fname xext]);
            delete([fpath '/m' fname xext]);
        catch ne_eo;
            warning( ...
                'neuroelf:DeleteFailed', ...
                'Error deleting IMG file: %s', ...
                ne_eo.message ...
            );
        end
    end
end
if ~keepnpar
    delete([fpath '/' fname '_seg_inv_sn.mat']);
    delete([fpath '/' fname '_seg_sn.mat']);
end

% change folder again
cd(opath);
