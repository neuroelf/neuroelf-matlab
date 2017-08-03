% FUNCTION ne_importvmr: convert an SPM image to a VMR
function varargout = ne_importvmr(varargin)

% Version:  v1.0
% Build:    16011511
% Date:     Jan-15 2016, 11:45 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2016, Jochen Weber
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

% block other callbacks (no graphics update, so this should work)
if ne_gcfg.c.incb
    return;
end
ne_gcfg.c.incb = true;

% request image file
[spmimgfile, spmimgpath] = uigetfile( ...
   {'*.hdr',  'Analyze 7.5 image files (*.hdr)'; ...
    '*.nii',  'NIftI image files (*.nii)'; ...
    '*.head', 'AFNI HEAD/BRIK files (*.head)'}, ...
    'Please select image file to convert to VMR...', ...
    '', 'MultiSelect', 'off');

% without valid name
if isequal(spmimgfile, 0) || ...
    isequal(spmimgpath, 0) || ...
    exist([spmimgpath '/' spmimgfile], 'file') ~= 2

    % reallow callbacks and return
    ne_gcfg.c.incb = false;
    return;
end

% patch path
spmimgfile = strrep(strrep([spmimgpath '/' spmimgfile], '\', '/'), '//', '/');

% get target filename
[vmrfile, vmrpath] = uiputfile({'*.vmr', 'BrainVoyager VMR files'}, ...
    'Please specify the target VMR filename...', '');

% and simply set to empty for bad path
if isequal(vmrfile, 0) || ...
    isequal(vmrpath, 0)
    vmrfile = '';
end

% patch target path as well
vmrfile = strrep(strrep([vmrpath '/' vmrfile], '\', '/'), '//', '/');

% try import
try
    tl{1} = [];
    th = importvmrfromanalyze(spmimgfile, 'cubic');
    if numel(th) ~= 1 || ...
       ~isxff(th)
        ne_gcfg.c.incb = false;
        return;
    end
    tl{1} = th;

    % save for non empty target name
    if ~isempty(vmrfile)
        th.SaveAs(vmrfile);
        disp('VMR successfully written:');
        disp(['  => ' th.FilenameOnDisk]);
    end

    % try to add to interface
    try
        ne_gcfg.c.incb = false;
        ne_openfile(0, 0, th, true);
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end

% on error
catch ne_eo;

    % store error
    ne_gcfg.c.lasterr = ne_eo;

    % clear any created object
    clearxffobjects(tl);

    % give an error
    errordlg('Error creating / saving VMR.', 'NeuroElf GUI - error', 'modal');

    % reallow callbacks and return
    ne_gcfg.c.incb = false;
    return;
end

% reallow callbacks
ne_gcfg.c.incb = false;
