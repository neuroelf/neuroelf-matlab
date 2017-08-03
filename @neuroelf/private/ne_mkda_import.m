% PUBLIC FUNCTION ne_mkda_import: import a PLP object for MKDA
function varargout = ne_mkda_import(varargin)

% Version:  v1.1
% Build:    16051820
% Date:     May-18 2016, 8:38 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, 2016, Jochen Weber
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
hFig = ne_gcfg.h.MKDA.MKDAFig;
ch = ne_gcfg.h.MKDA.h;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% blocked?
if any(strcmp(ne_gcfg.c.blockcb, 'mkda_load'))
    return;
end
ne_gcfg.c.blockcb{end+1} = 'mkda_load';

% get content of PLP dropdown
plps = ch.PLPs;
plpud = plps.UserData;
plpst = plps.String;
if ~iscell(plpst)
    plpst = cellstr(plpst);
end

% request filename for import
mfp = hFig.Pointer;
try
    [plpifile, plpipath] = uigetfile({'*.mat;*.xls;*.xlsx;*.csv;*.txt', ...
        'MKDA database files (*.mat, *.xls, *.xlsx, *.csv, *.txt)'}, ...
        'Please select a file to import...');
    if isequal(plpifile, 0) || ...
        isequal(plpipath, 0) || ...
        isempty(plpifile)
        ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'mkda_load')) = [];
        return;
    end
    if isempty(plpipath)
        plpipath = pwd;
    end
    plpffile = strrep(strrep([plpipath '/' plpifile], '\', '/'), '//', '/');
    hFig.Pointer = 'watch';
    drawnow;
    plp = importplp(plpffile);

% give warning
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    uiwait(warndlg('Error loading PLP file.', 'NeuroElf GUI - error', 'modal'));
    hFig.Pointer = mfp;
end

% add saved flag
plp.RunTimeVars.Saved = false;

% add to lists
plpud(end+1, :) = {cell(0, 2), [], plp};
if numel(plpst) == 1 && ...
    strcmpi(plpst{1}, '<no plp loaded>')
    plpst(1) = [];
end
if numel(plpifile) <= 36
    plpst{end+1} = plpifile;
else
    plpst{end+1} = [plpifile(1:14) '...' plpifile(end-17:end)];
end
plps.String = plpst(:);
plps.UserData = plpud;
plps.Value = numel(plpst);

% unblock
ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'mkda_load')) = [];

% enable PLP menu
hFig.SetGroupEnabled('PLPOK', 'on');

% and make current PLP
ne_mkda_setplp;

% revert pointer
hFig.Pointer = mfp;
