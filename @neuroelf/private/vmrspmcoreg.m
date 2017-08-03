function vmrtrf = vmrspmcoreg(srcvmr, refvmr)
% vmrspmcoreg  - apply SPM's coregistration to two VMRs
%
% FORMAT:       vmrtrf = vmrspmcoreg(srcvmr, refvmr)
%
% Input fields:
%
%       srcvmr      VMR to coregister (primary VMR)
%       refvmr      reference (stationary/secondary VMR)
%
% Output fields:
%
%       vmrtrf      BV compatible TRF object

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:31 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2014, 2016, Jochen Weber
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

% is SPM available
try
    sver = lower(spm('ver'));
    if ~any(strcmp(sver, {'spm5', 'spm8', 'spm12b', 'spm12'}))
        error( ...
            'neuroelf:SPMError', ...
            'Required SPM version not on the path.' ...
        );
    end
    sver = str2double(regexprep(sver, '^spm(\d+).*$', '$1'));
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    error( ...
        'neuroelf:SPMError', ...
        'SPM not found on the path.' ...
    );
end

% load appropriate template
try
    if sver == 5
        jobs = neuroelf_file('p', 'spm5_coreg');
    else
        jobs = neuroelf_file('p', 'spm8_coreg');
    end
    jobs = jobs.jobs;
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    error( ...
        'neuroelf:PackageError', ...
        'Required coregistration job not found.' ...
    );
end

% argument check (and/or selection)
clearo = cell(1, 2);
if nargin < 1 || ...
   ~isxff(srcvmr, 'vmr')
    try
        srcvmr = xff('*.vmr', 'Please select the source (primary) VMR...');
        clearo{1} = srcvmr;
        if ~isxff(srcvmr, 'vmr')
            error('NOT_A_VMR');
        end
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        clearxffobjects(clearo);
        error( ...
            'neuroelf:BadArgument', ...
            'Invalid VMR selected.' ...
        );
    end
end
if nargin < 2 || ...
   ~isxff(refvmr, 'vmr')
    try
        refvmr = xff('*.vmr', 'Please select the reference (secondary) VMR...');
        clearo{2} = refvmr;
        if ~isxff(refvmr, 'vmr')
            error('NOT_A_VMR');
        end
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        clearxffobjects(clearo);
        error( ...
            'neuroelf:BadArgument', ...
            'Invalid VMR selected.' ...
        );
    end
end

% progress in console
disp('Writing VMRs -> Analyze...');

% write both as temporary files
srcimg = [tempname(tempdir) '.img'];
refimg = [tempname(tempdir) '.img'];
srcvmr.WriteAnalyzeVol(srcimg);
refvmr.WriteAnalyzeVol(refimg);

% keep filenames for TRF
srcfile = srcvmr.FilenameOnDisk;
if isempty(srcfile)
    srcfile = '<none>';
end
reffile = refvmr.FilenameOnDisk;
if isempty(reffile)
    reffile = '<none>';
end

% clear temporary object
clearxffobjects(clearo);

% put files into structure
if sver == 5
    jobs{1}.spatial{1}.coreg{1}.estimate.source{1} = srcimg;
    jobs{1}.spatial{1}.coreg{1}.estimate.ref{1} = refimg;
    spm_defaults;
else
    jobs{1}.spm.spatial.coreg.estimate.source{1} = srcimg;
    jobs{1}.spm.spatial.coreg.estimate.ref{1} = refimg;
    spm('defaults', 'FMRI');
end

% try to run job
disp('Running SPM coregistration...');
spm_jobman('run', jobs);

% read both headers
clearo = cell(1, 2);
try
    srchdr = xff(strrep(srcimg, '.img', '.hdr'));
    clearo{1} = srchdr;
    refhdr = xff(strrep(refimg, '.img', '.hdr'));
    clearo{2} = refhdr;
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    clearxffobjects(clearo);
    error( ...
        'neuroelf:FileError', ...
        'Error opening header files after coregistration.' ...
    );
end

% get required transform (do the +0.5/-0.5 thing to cope with different
% origin interpretations!!)
trf = spmtrf([-0.5, -0.5, -0.5]) * ...
    refhdr.CoordinateFrame.Trf * inv(srchdr.CoordinateFrame.Trf) * ...
    spmtrf([0.5, 0.5, 0.5]);

% clear objects
clearxffobjects(clearo);

% delete files
delete(srcimg);
delete(refimg);
delete(strrep(srcimg, '.img', '.hdr'));
delete(strrep(refimg, '.img', '.hdr'));

% create new TRF
disp('Creating TRF...');
vmrtrf = xff('new:trf');

% make settings
vmrtrf.TransformationType = 2;
vmrtrf.CoordinateSystem = 0;
vmrtrf.ExtraVMRTransf = 0;
vmrtrf.ExtraVMRTransfValues = eye(4);
vmrtrf.SourceFile = srcfile;
vmrtrf.TargetFile = reffile;

% make sure TRF.TFMatrix is correct !
trf = [trf([2, 3, 1], [2, 3, 1, 4]); 0, 0, 0, 1];
trf(1:3, 4) = -trf(1:3, 4);
vmrtrf.TFMatrix = trf;
