% PUBLIC FUNCTION ne_exportnii: export the selected slice/stats file
function varargout = ne_exportnii(varargin)

% Version:  v0.9c
% Build:    11050712
% Date:     Apr-24 2011, 9:40 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, Jochen Weber
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
ch = ne_gcfg.h;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% blocked
if ne_gcfg.c.incb || ...
    any(strcmp(ne_gcfg.c.blockcb, 'exportnii'))
    return;
end

% inputs?
if nargin < 3 || ...
   ((~ischar(varargin{3}) || ...
     ~any(strcmpi(varargin{3}(:)', {'slicevar', 'statsvar'}))) && ...
    (numel(varargin{3}) ~= 1 || ...
     ~isxff(varargin{3}, ...
        {'ava', 'cmp', 'dmr', 'fmr', 'glm', 'head', 'msk', 'vmp', 'vmr', 'vtc'})))
    return;
end
if ischar(varargin{3})
    switch (lower(varargin{3}(:)'))
        case {'slicevar'}
            evar = cc.SliceVar;
        case {'statsvar'}
            evar = cc.StatsVar;
    end
    if ~isxff(evar, {'ava', 'cmp', 'dmr', 'fmr', 'glm', 'head', 'msk', 'vmp', 'vmr', 'vtc'})
        return;
    end
else
    evar = varargin{3};
end

% get filename
efile = evar.FilenameOnDisk;
if isempty(efile)
    efile = 'exported.nii';
end
[efld, efile] = fileparts(efile);
efile = [efile '.nii'];
[efile, efld] = uiputfile( ...
    {'*.nii', 'NIftI file (*.nii)'; ...
     '*.hdr', 'Analyze 7.5 file (*.hdr)'}, ...
    'Export object to NII...', efile);
if isequal(efile, 0) || ...
    isequal(efld, 0)
    return;
end

% block callbacks, pointer
ne_gcfg.c.incb = true;
ne_gcfg.c.blockcb{end+1} = 'exportnii';
mfp = ch.MainFig.Pointer;
ch.MainFig.Pointer = 'watch';
drawnow;

% get var to save
try
    evar.ExportNifti([efld '/' efile]);
catch ne_eo;
    uiwait(warndlg(['Error exporting object to NII: ' ne_eo.message], ...
        'NeuroElf - warning', 'modal'));
end

% unset pointer, block, etc.
ch.MainFig.Pointer = mfp;
drawnow;
ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'exportnii')) = [];
ne_gcfg.c.incb = false;
