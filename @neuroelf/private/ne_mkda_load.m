% PUBLIC FUNCTION ne_mkda_load: load a PLP object for MKDA
function varargout = ne_mkda_load(varargin)

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:38 AM EST
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

% request plp
mfp = hFig.Pointer;
hFig.Pointer = 'watch';
drawnow;
plp = [];
try
    if nargin < 3 || ...
       ~ischar(varargin{3}) || ...
        isempty(varargin{3})
        if nargin > 2 && ...
            numel(varargin{3}) == 1 && ...
            isxff(varargin{3}, 'plp')
            plp = varargin{3};
        else
            plp = xff('*.plp');
        end
    else
        plp = xff(varargin{3}(:)', 'plp');
    end
    if isempty(plp)
        ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'mkda_load')) = [];
        hFig.Pointer = mfp;
        return;
    end
    if numel(plp) ~= 1 || ...
       ~isxff(plp, 'plp')
        error( ...
            'neuroelf:BadXFFOutput', ...
            'Not a valid PLP object.' ...
        );
    end

% give warning
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    if numel(plp) == 1 && ...
        isxff(plp, true)
        plp.ClearObject;
    end
    uiwait(warndlg('Error loading PLP file.', 'NeuroElf GUI - error', 'modal'));
    ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'mkda_load')) = [];
    hFig.Pointer = mfp;
    return;
end

% add saved flag
plp.RunTimeVars.Saved = true;

% add to lists
plpud(end+1, :) = {cell(0, 2), [], plp};
[plpp, plpf, plpe] = fileparts(plp.FilenameOnDisk);
if numel(plpst) == 1 && ...
    strcmpi(plpst{1}, '<no plp loaded>')
    plpst(1) = [];
end
if numel(plpf) <= 32
    plpst{end+1} = [plpf plpe];
else
    plpst{end+1} = [plpf(1:14) '...' plpf(end-13:end) plpe];
end
plps.String = plpst(:);
plps.UserData = plpud;
plps.Value = numel(plpst);

% unblock
ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'mkda_load')) = [];

% and make current PLP
ne_mkda_setplp;

% revert pointer
hFig.Pointer = mfp;
