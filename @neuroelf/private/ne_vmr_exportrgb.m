% FUNCTION ne_vmr_exportrgb: export RGB NIftI image
function varargout = ne_vmr_exportrgb(varargin)

% Version:  v0.9d
% Build:    14071116
% Date:     Jul-11 2014, 4:40 PM EST
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

% global variable
global ne_gcfg;
cc = ne_gcfg.fcfg;

% initialize output
varargout = cell(1, nargout);

% test SliceVar
svar = cc.SliceVar;
if numel(svar) ~= 1 || ...
   ~isxff(svar, 'vmr')
    return;
end

% get target filename
vmrfile = svar.FilenameOnDisk;
if isempty(vmrfile)
    vmrfile = [pwd '/unsaved.vmr'];
end
[vmrpath, vmrfile] = fileparts(vmrfile);
[niifile, niipath] = uiputfile({'*.nii', 'NIftI files (*.nii)'}, ...
    'Please choose a target file location.', [vmrfile '.nii']);
if isequal(niifile, 0) || ...
    isequal(niipath, 0) || ...
    isempty(niifile)
    return;
end
if isempty(niipath)
    niipath = vmrpath;
end

% export
try
    niifile = [niipath '/' niifile];
    svar.ExportNIftI(niifile, true);
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    return;
end

% read the file
svar = ne_openfile(0, 0, niifile, true);

% adapt content
svar.VoxelDataRGBA = repmat(svar.VoxelData(:, :, :), [1, 1, 1, 1, 3]);
svar.VoxelData = [];
svar.ImgDim.DataType = 128;
svar.Save;
if isfield(svar.RunTimeVars, 'UndoBuffer')
    svar.RunTimeVars = rmfield(svar.RunTimeVars, 'UndoBuffer');
    svar.SaveRunTimeVars;
end

% update UI
ne_setslicepos;
drawnow;
