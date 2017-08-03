% FUNCTION ne_renamedicom: bring up rename dicom dialog
function varargout = ne_renamedicom(varargin)

% Version:  v1.0
% Build:    16010821
% Date:     Jan-08 2016, 9:29 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2014, 2016, Jochen Weber
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

% only allow one instance at a time
if any(strcmp('renamedicom', ne_gcfg.c.blockcb))
    return;
end
ne_gcfg.c.blockcb{end+1} = 'renamedicom';

% hour glass
mf = ne_gcfg.h.MainFig;
mfp = mf.Pointer;
mf.Pointer = 'watch';
cprog = ne_progress(0, 0, {true, 0, 'renamedicom'});
prb = ne_gcfg.h.Progress;

% call renamedicom with the progress bar handle
try
    p = renamedicom(struct('pbar', prb));
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    uiwait(warndlg(ne_eo.message, 'NeuroElf - error', 'modal'));
    mf.Pointer = mfp;
    ne_progress(0, 0, cprog);
    ne_gcfg.c.blockcb(strcmp('renamedicom', ne_gcfg.c.blockcb)) = [];
    return;
end

% packed renaming?
if ~isempty(p)

    % all new filenames are "local" filenames
    lf = all(cellfun('isempty', regexp(p(:, 2), '[\\\/]')));
    if lf

        % ask for what else to do
        action = questdlg('New filenames are relative. Where to store files?', ...
            'NeuroElf - user input', 'Original folder(s)', 'Current folder', ...
            'Specify folder', 'Specify folder');
        if ~ischar(action) || ...
            isempty(action)
            action = 's';
        else
            action = lower(action(1));
        end
    else
        action = 'o';
    end

    % specify target folder
    if action == 's'
        tf = uigetdir(pwd, 'Please select folder to store DICOM files...');
        if ~ischar(tf) || ...
            isempty(tf)
            mf.Pointer = mfp;
            ne_progress(0, 0, cprog);
            ne_gcfg.c.blockcb(strcmp('renamedicom', ne_gcfg.c.blockcb)) = [];
            return;
        end
    elseif action == 'c'
        tf = pwd;
    end

    % write out files
    fc = 1;
    nc = size(p, 1);
    lpt = now;
    lpd = 1 / 86400;
    prb.Progress(0, sprintf('Writing out %d DICOM files...', nc));
    try
        if action == 'o'
            for fc = 1:nc
                tf = fileparts(p{fc, 1});
                [nullf, tfile, text] = fileparts(p{fc, 2});
                if (now - lpt) >= lpd
                    lpt = now;
                    prb.Progress((fc - 1) / nc, sprintf(...
                        'Writing out ''%s'' (%d/%d)...', [tfile text], fc, nc));
                end
                binwrite([tf '/' tfile text], p{fc, 3});
            end
        else
            for fc = 1:nc
                [nullf, tfile, text] = fileparts(p{fc, 2});
                if (now - lpt) >= lpd
                    lpt = now;
                    prb.Progress((fc - 1) / nc, sprintf(...
                        'Writing out ''%s'' (%d/%d)...', [tfile text], fc, nc));
                end
                binwrite([tf '/' tfile text], p{fc, 3});
            end
        end
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
        uiwait(warndlg(sprintf('Error writing out file %d/%d: %s.', ...
            fc, nc, ne_eo.message), 'NeuroElf - error', 'modal'));
        mf.Pointer = mfp;
        ne_progress(0, 0, cprog);
        ne_gcfg.c.blockcb(strcmp('renamedicom', ne_gcfg.c.blockcb)) = [];
        return;
    end
end

% pointer and progress bar reset
mf.Pointer = mfp;
ne_progress(0, 0, cprog);

% reset the callback blocker
ne_gcfg.c.blockcb(strcmp('renamedicom', ne_gcfg.c.blockcb)) = [];
