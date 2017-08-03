% FUNCTION ne_vmp_writeana: write Analyze version of selected volumes
function varargout = ne_vmp_writeana(varargin)

% Version:  v0.9b
% Build:    11050712
% Date:     Apr-10 2011, 4:52 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, Jochen Weber
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

% only valid if single VMP and at least one SRF in scenery is selected
vmp = ne_gcfg.fcfg.StatsVar;
if numel(vmp) ~= 1 || ...
   ~isxff(vmp, 'vmp')
    return;
end

% get maps to create MSK from
mapsel = ch.StatsVarMaps.Value;
if isempty(mapsel)
    return;
end

% get the filename pattern
fp = ne_gcfg.c.ini.VMP.WriteAnaPattern;

% get filename of VMP
[vpath, vfile] = fileparts(vmp.FilenameOnDisk);
if isempty(vpath)
    vpath = pwd;
end

% replace common pattern ?
if ~isempty(strfind(fp, '$f'))

    % check name is valid
    if isempty(vfile)
        vfile = 'unnamed';
    end

    % replace
    fp = strrep(fp, '$f', vfile);
end

% set pointer
pt = ch.MainFig.Pointer;
ch.MainFig.Pointer = 'watch';
drawnow;

% with error handling...
try
    % iterate over maps
    for mc = mapsel(:)'

        % need to replace map name ?
        if ~isempty(strfind(fp, '$m'))

            % get map name as label
            mn = makelabel(vmp.Map(mc).Name);

            % and replace
            vfile = strrep(fp, '$m', mn);

        % otherwise
        else

            % simply copy name
            vfile = fp;
        end

        % map number?
        if ~isempty(regexpi(vfile, '\%\d+d'))

            % replace
            vfile = sprintf(vfile, mc);
        end

        % save as
        vmp.WriteAnalyzeVol(mc, sprintf('%s/%s', vpath, vfile));
    end

% error?
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;

    uiwait(warndlg(sprintf('Error writing Analyze file for map %d', mc), ...
        'NeuroElf - error', 'modal'));
end

% reset pointer
ch.MainFig.Pointer = pt;
drawnow;
