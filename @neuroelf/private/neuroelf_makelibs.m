function neuroelf_makelibs(hPrg, hBar, pRng)
% neuroelf_makelibs  - compile MEX functions (mex/compiler needed!)
%
% FORMAT:       neuroelf_makelibs;
%
% No input/output fields.
%
% Note: this should create the following (system dependent) files:
%
%       applyspmsnc.EXT                 - apply SPM normalization
%       clustercoordsc.EXT              - cluster a volume
%       condsum.EXT                     - conditional sum of values
%       conv3d.EXT                      - 3D convolution
%       cov_nd.EXT                      - covariance along last dim
%       findfirst.EXT                   - find first (last) true value
%       flexinterpn.EXT                 - flexible up-to 4D interpolation
%       flexmask.EXT                    - flexible masking
%       floodfill3c.EXT                 - flood-filling of slice/volume
%       gluetostringc.EXT               - concatenate strings
%       histcount.EXT                   - histogram counting (2D/weighted)
%       hsvconv.EXT                     - HSV/RGB colorspace conversion
%       indexarray.EXT                  - array indexing with ranges
%       indexarraynb.EXT                - out-of-bounds-ok array indexing
%       isinfnan.EXT                    - test single/double for Inf/NaN's
%       joinlayers.EXT                  - join layers (transimg)
%       kendtaupairsign.EXT             - compute Kendall Tau pair sign
%       limitrangec.EXT                 - limit numeric data to range
%       mesh_dist.EXT                   - compute mesh point-pair distances
%       mesh_morph.EXT                  - mesh morphing
%       mesh_neighborsarray.EXT         - generate neighbors array
%       mesh_normals.EXT                - mesh normals computation
%       mesh_reconstruct.EXT            - marching-cube reconstruction
%       mesh_trianglestoneighbors.EXT   - triangle-to-neighbors conversion
%       mesh_trimapmesh.EXT             - mesh mapping vertex-to-triangle
%       minmaxmean.EXT                  - find min, max, compute mean/std
%       normcdfc.EXT                    - normal CDF C-implementation
%       psetdists.EXT                   - compute all pair-wise distances
%       renderlayers.EXT                - render layers (transimg)
%       renderv3dub.EXT                 - update-buffer MEX for renderv3d
%       renderv3dxia.EXT                - access indices MEX for renderv3d
%       setsparseval.EXT                - set values of sparse array
%       splittocellc.EXT                - split string into cell array
%       sprintfbx.EXT                   - print double/singles as hex (%bx)
%       svmpredictc2x.EXT               - libSVM v2.89 C-implementation
%       svmpredictc3x.EXT               - libSVM v3.18 C-implementation
%       svmreadc3x.EXT                  - libSVM v3.18 C-implementation
%       svmtrainc2x.EXT                 - libSVM v2.89 C-implementation
%       svmtrainc3x.EXT                 - libSVM v3.18 C-implementation
%       svmwritec3x.EXT                 - libSVM v3.18 C-implementation
%       threshlutc.EXT                  - threshold values LUT -> RGB
%       threshmapc.EXT                  - threshold map values
%       transmul.EXT                    - transposed matrix multiplication
%       u8str2double.EXT                - convert 8-bit strings to double
%       varc.EXT                        - var C-implementation
%       xffsrfparseneighborsc.EXT       - xff SRF neighbors parsing helper
%       xffsrfwriteneighborsc.EXT       - xff SRF neighbors writing helper
%
% If you wish to help the distributors of NeuroElf, please compile
% these files with Matlab v7.0 (R14) on your desired platform and
% send in the binary results to
%
% info@neuroelf.net
%
% Thank you!

% Version:  v1.0
% Build:    14121416
% Date:     Dec-14 2014, 4:43 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 -2014, Jochen Weber
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

% arguments
if nargin < 1
    hPrg = [];
end
if nargin < 2
    hBar = [];
end
if nargin < 3
    pRng = [0, 1];
end
pMin = pRng(1);
pDif = pRng(2) - pMin + eps;

% get old and new pwd
opwd = pwd;
npwd = [neuroelf_path filesep '@neuroelf' filesep 'private'];

% change to path
eoc = 0;
ismatlab = (~isempty(ver('matlab')) || isempty(ver('octave')));
try
    cd(npwd);
    cfiles = findfiles(npwd, '*.c', 'depth=1');
    ler = cell(numel(cfiles), 1);
    mext = mexext;
    for cc = 1:numel(cfiles)
        [mfp, mfn, mfe] = fileparts(cfiles{cc});
        if ~isempty(regexpi(mfn, '^svm'))
            continue;
        end
        try
            mexx = dir([mfn '.' mext]);
            if numel(mexx) == 1 && ...
                mexx.isdir == 0 && ...
                mexx.bytes >= 4096
                instprog(hPrg, hBar, sprintf('MEX found: %s', mfn), ...
                    pMin + pDif * (cc / (numel(cfiles) + 4)));
                continue;
            end
            instprog(hPrg, hBar, sprintf('Compiling %s...', cfiles{cc}), ...
                pMin + pDif * (cc / (numel(cfiles) + 4)), ...
                ['Compiling MEX files: ' mfn '...']);
            if ismatlab
                if ~any(strcmpi(mfn, {'setsparseval'}))
                    mex('-O', [mfn mfe]);
                else
                    mex('-O', '-largeArrayDims', [mfn mfe]);
                end
            else
                mex([mfn mfe]);
                try
                    delete([mfn '.o']);
                catch ne_eo;
                    neuroelf_lasterr(ne_eo);
                end
            end
            if exist([mfn '.' mext], 'file') < 1
                eoc = eoc + 1;
                break;
            end
        catch ne_eo;
            eoc = eoc + 1;
            ler{eoc} = ne_eo.message;
        end
    end

    % libSVM files need special treatment
    cd(npwd);
    is64bit = ~isempty(regexpi(mext, '64'));
    instprog(hPrg, hBar, ...
        'Compiling libSVM v2.x/3.x files...', ...
        pMin + pDif * (numel(cfiles) / (numel(cfiles) + 3)));
    if exist([pwd filesep 'svmpredictc2x.' mext], 'file') > 0
        instprog(hPrg, hBar, 'MEX found: svmpredictc2x');
    else
        instprog(hPrg, hBar, ...
            ['Compiling MEX files: ' pwd filesep 'svmpredictc2x...']);
        if is64bit
            mex('-O', '-largeArrayDims', 'svmpredictc2x.c', 'svm2x.cpp');
        elseif ismatlab
            mex('-O', 'svmpredictc2x.c', 'svm2x.cpp');
        else
            mex('svmpredictc2x.c', 'svm2x.cpp');
            try
                delete('svmpredictc2x.o');
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
            end
            try
                delete('svm2x.o');
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
            end
        end
    end
    if exist([pwd filesep 'svmtrainc2x.' mext], 'file') > 0
        instprog(hPrg, hBar, 'MEX found: svmtrainc2x');
    else
        instprog(hPrg, hBar, ...
            ['Compiling MEX files: ' pwd filesep 'svmtrainc2x...']);
        if is64bit
            mex('-O', '-largeArrayDims', 'svmtrainc2x.c', 'svm2x.cpp');
        elseif ismatlab
            mex('-O', 'svmtrainc2x.c', 'svm2x.cpp');
        else
            mex('svmtrainc2x.c', 'svm2x.cpp');
            try
                delete('svmtrainc2x.o');
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
            end
            try
                delete('svm2x.o');
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
            end
        end
    end
    if exist([pwd filesep 'svmpredictc3x.' mext], 'file') > 0
        instprog(hPrg, hBar, 'MEX found: svmpredictc3x');
    else
        instprog(hPrg, hBar, ...
            ['Compiling MEX files: ' pwd filesep 'svmpredictc3x...']);
        if is64bit
            mex('-O', '-largeArrayDims', 'svmpredictc3x.c', 'svm3x.cpp', 'svm3xmodel.c');
        elseif ismatlab
            mex('-O', 'svmpredictc3x.c', 'svm3x.cpp', 'svm3xmodel.c');
        else
            mex('svmpredictc3x.c', 'svm3x.cpp', 'svm3xmodel.c');
            try
                delete('svmpredictc3x.o');
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
            end
            try
                delete('svm3x.o');
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
            end
            try
                delete('svm3xmodel.o');
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
            end
        end
    end
    if exist([pwd filesep 'svmreadc3x.' mext], 'file') > 0
        instprog(hPrg, hBar, 'MEX found: svmreadc3x');
    else
        instprog(hPrg, hBar, ...
            ['Compiling MEX files: ' pwd filesep 'svmreadc3x...']);
        if is64bit
            mex('-O', '-largeArrayDims', 'svmreadc3x.c');
        elseif ismatlab
            mex('-O', 'svmreadc3x.c');
        else
            mex('svmreadc3x.c');
            try
                delete('svmreadc3x.o');
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
            end
        end
    end
    if exist([pwd filesep 'svmtrainc3x.' mext], 'file') > 0
        instprog(hPrg, hBar, 'MEX found: svmtrainc3x');
    else
        instprog(hPrg, hBar, ...
            ['Compiling MEX files: ' pwd filesep 'svmtrainc3x...']);
        if is64bit
            mex('-O', '-largeArrayDims', 'svmtrainc3x.c', 'svm3x.cpp', 'svm3xmodel.c');
        elseif ismatlab
            mex('-O', 'svmtrainc3x.c', 'svm3x.cpp', 'svm3xmodel.c');
        else
            mex('svmtrainc3x.c', 'svm3x.cpp', 'svm3xmodel.c');
            try
                delete('svmtrainc3x.o');
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
            end
            try
                delete('svm3x.o');
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
            end
            try
                delete('svm3xmodel.o');
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
            end
        end
    end
    if exist([pwd filesep 'svmwritec3x.' mext], 'file') > 0
        instprog(hPrg, hBar, 'MEX found: svmwritec3x');
    else
        instprog(hPrg, hBar, ...
            ['Compiling MEX files: ' pwd filesep 'svmwritec3x...']);
        if is64bit
            mex('-O', '-largeArrayDims', 'svmwritec3x.c');
        elseif ismatlab
            mex('-O', 'svmwritec3x.c');
        else
            mex('svmwritec3x.c');
            try
                delete('svmwritec3x.o');
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
            end
        end
    end
catch ne_eo;
    cd(opwd);
    error( ...
        'neuroelf:MEXError', ...
        'Custom error compiling MEX functions: %s', ...
        ne_eo.message ...
    );
end

% change pwd back
cd(opwd);

% specific error message
if eoc > 0
    error( ...
        'neuroelf:MEXError', ...
        'Error compiling at least one function.' ...
    );
end


% sub function for progress
function instprog(hedt, hbar, txt, p, ptxt)
if numel(hedt) == 1
    if ~isempty(txt)
        str = hedt.String;
        if ischar(str) && ...
            size(str, 1) > 1
            str = cellstr(str);
        elseif ischar(str)
            str = splittocell(str, char(10));
        end
        str = [str(:); lsqueeze(splittocell(txt, char(10)))];
        hedt.String = str;
        hedt.Value = numel(str);
        hedt.ListboxTop = max(1, numel(str) - 15);
    end
    if nargin > 4
        hbar.Progress(p, ptxt);
    elseif nargin > 3
        hbar.Progress(p);
    else
        drawnow;
    end
elseif isempty(hedt)
    disp(txt);
end
