% PUBLIC FUNCTION ne_mkda_save: save a PLP object from MKDA
function varargout = ne_mkda_save(varargin)

% Version:  v0.9c
% Build:    11120709
% Date:     Nov-08 2011, 3:23 PM EST
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
hFig = ne_gcfg.h.MKDA.MKDAFig;
ch = ne_gcfg.h.MKDA.h;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% saveas
if nargin < 3 || ...
   ~ischar(varargin{3}) || ...
   ~strcmpi(varargin{3}(:)', 'saveas')
    saveas = false;
else
    saveas = true;
end

% for new files, only allow saveas
plp = ne_gcfg.fcfg.plp;
if numel(plp) ~= 1 || ...
   ~isxff(plp, 'plp')
    return;
end
if isempty(plp.FilenameOnDisk)
    saveas = true;
end

% depending on saveas flag
mfp = hFig.Pointer;
try
    plp.RunTimeVars.AutoSave = true;
    hFig.Pointer = 'watch';
    drawnow;
    if saveas
        plp.SaveAs;
        fnames = ch.PLPs.String;
        if ~iscell(fnames)
            fnames = cellstr(fnames);
        end
        [plpp, plpf, plpe] = fileparts(plp.FilenameOnDisk);
        if numel(plpf) <= 32
            fnames{ch.PLPs.Value} = [plpf plpe];
        else
            fnames{ch.PLPs.Value} = [plpf(1:14) '...' plpf(end-13:end) plpe];
        end
        ch.PLPs.String = fnames;
    else
        plp.Save;
    end
    plp.RunTimeVars.Saved = true;
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    uiwait(warndlg('Error saving PLP file.', 'NeuroElf GUI - error', 'modal'));
end
hFig.Pointer = mfp;
