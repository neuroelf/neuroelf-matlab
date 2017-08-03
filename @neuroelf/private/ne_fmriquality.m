% FUNCTION ne_fmriquality: invoke fmriquality after requesting some files
function varargout = ne_fmriquality(varargin)

% Version:  v1.0
% Build:    16010821
% Date:     Jan-08 2016, 9:29 PM EST
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
ch = ne_gcfg.h;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% only allow one instance at a time
if any(strcmp('fmriquality', ne_gcfg.c.blockcb))
    return;
end
ne_gcfg.c.blockcb{end+1} = 'fmriquality';

% if needed, request filename(s) to work on
if nargin < 3 || ...
   ~iscell(varargin{3}) || ...
    isempty(varargin{3}) || ...
   ~all(cellfun(@ischar, varargin{3}(:)))

    [fmrifiles, fmrifpath, fmrifi] = uigetfile( ...
        {'*.hdr;*.nii', 'Analyze images (*.hdr, *.nii)'; ...
         '*.fmr',       'Functional MR files (*.fmr)'; ...
         '*.vtc',       'Volume Time-Course files (*.vtc)'; ...
         '*.fq',        'fMRI Quality files (*.fq)'}, ...
        'Please select the file(s) of one fMRI run...', 'MultiSelect', 'on');
    if isequal(fmrifiles, 0) || ...
        isequal(fmrifpath, 0) || ...
        isempty(fmrifiles)
        ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'fmriquality')) = [];
        return;
    end
    if ~iscell(fmrifiles)
        fmrifiles = cellstr(fmrifiles);
    end
    for fc = 1:numel(fmrifiles)
        fmrifiles{fc} = [fmrifpath fmrifiles{fc}];
    end
    if fmrifi == 4
        for fc = 1:numel(fmrifiles)
            try
                q = load(fmrifiles{fc}, '-mat');
                ne_gcfg.h.fChild(end+1) = fmriqasheet(q.fq);
            catch ne_eo;
                 ne_gcfg.c.lasterr = ne_eo;
            end
        end
        ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'fmriquality')) = [];
        return;
    end
else
    fmrifiles = varargin{3}(:);
end

% load first file to determine type
try
    fmrifirst = fmrifiles{1};
    fmrifiles{1} = xff(fmrifirst);
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    uiwait(warndlg('Error reading first file for fmriquality.', ...
        'NeuroElf GUI - Info', 'modal'));
    ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'fmriquality')) = [];
    return;
end

% for FMR and VTCs, reject rest (if any)
if any(strcmpi(fmrifiles{1}.Filetype, {'fmr', 'vtc'}))
    if fmrifi > 1 && ...
        numel(fmrifiles) > 1
        uiwait(warndlg('Only first FMR/VTC file considered!', 'NeuroElf GUI - Info', 'modal'));
        fmrifiles(2:end) = [];
    end
    fmrifiles = fmrifiles{1};

% for Analyze
else

    % put back the filename
    fmrifiles{1}.ClearObject;
    fmrifiles{1} = fmrifirst;
end

% give feedback
if ne_gcfg.c.echo
    fmrioa = 'struct(''motcor'', true, ''qasheet'', false)';
    if numel(fmrifiles) < 4
        fmrifa = any2ascii(fmrifiles);
        ne_echo({'qasheet = fmriquality(%s, %s)', fmrifa, fmrioa});
    else
        ne_echo({'qasheet = fmriquality(fmrifiles, %s);', fmrioa});
    end
end

% set progress bar visible and call fmriquality
cprog = ne_progress(0, 0, {true, 0, 'fmriquality'});
try
    q = fmriquality(fmrifiles, struct( ...
        'pbar', ch.Progress, 'motcor', true, 'qasheet', false));
    save([fmrifirst(1:end-4) '_fmriquality.mat'], 'q');
    ne_gcfg.h.fChild(end+1) = fmriqasheet(q);
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    warning( ...
        'neuroelf:FileError', ...
        'Error saving fmriquality output to MAT file.' ...
    );
end
ne_progress(0, 0, cprog);

% clear object?
if ~iscell(fmrifiles)
    fmrifiles.ClearObject;
end

% unblock next call
ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'fmriquality')) = [];
