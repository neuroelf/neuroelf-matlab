function varargout = ne_closesatwindow(varargin)
% ne_closesatwindow  - close satellite window
%
% FORMAT:       ne_closesatwindow(SRC, EVT, windowid)
%
% Input fields:
%
%       SRC, EVT    Matlab handle callback inputs (discarded)
%       windowid    1x8 string identifying window to be closed
%
% No output fields.
%
% Example using scripting:
%
% % show surface window
% neuroelf_gui('showpage', 3);
%
% % iterate over time (morphing) values
% for time = 0:0.01:1
%
%     % set position
%     neuroelf_gui('setsurfpos', '', {150, 15, [0, -30, 0], 1.25, time});
%
%     % undock
%     [sat, ~, satid] = neuroelf_gui('undock');
%
%     % set correct aspect ratio, etc.
%     neuroelf_gui('satresize', satid, [720, 480]);
%
%     % screenshot
%     neuroelf_gui('screenshot', satid, ...
%         sprintf('scene_time%03d.png', 100 * time), 'high-q');
%
%     % **CLOSE WINDOW**
%     neuroelf_gui('closesatwindow', satid);
% end

% Version:  v1.1
% Build:    16052610
% Date:     May-26 2016, 10:34 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2016, Jochen Weber
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

% check input
if nargin < 3 || ~ischar(varargin{3}) || ~isfield(ne_gcfg.cc, varargin{3})

    % presumably correct input
    if nargin > 2 && ischar(varargin{3}) && numel(varargin{3}) == 8

        % try to assess handle
        try
            figh = gcbf;
            figp = get(figh);
            if ~strcmpi(figp.Type, 'figure') || isempty(regexpi(figp.Name, '^neuroelf'))
                return;
            end

            % presumably, we lost contact (reloaded figure?)

            % unset DeleteFcn
            set(figh, 'DeleteFcn', '');
            set(get(figh, 'Children'), 'DeleteFcn', '');

            % then delete
            delete(figh);
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end
    end
    return;
end
ch = ne_gcfg.cc.(varargin{3});

% PLP object
if isfield(ch.Config, 'plp') && numel(ch.Config.plp) == 1 && isxff(ch.Config.plp, 'plp')
    ch.Config.plp.ClearObject;
    ne_gcfg.cc.(varargin{3}).Config.plp = [];
end

% delete satellite window and remove from children
try
    lastpos = ch.Satellite.Position;
    ch.Satellite.Delete;

    % delete transimg objects?
    if strcmp(ch.Config.sattype, 'slice')
        delete(ch.Config.tio.imSag);
        delete(ch.Config.tio.imCor);
        delete(ch.Config.tio.imTra);
        delete(ch.Config.tio.imSlZ);

        % and keep track of position
        ne_gcfg.c.ini.Satellites.Position = lastpos(1:2);

    % for plot window
    elseif strcmp(ch.Config.sattype, 'betaplot')

        % record position
        ne_gcfg.c.ini.Children.BetaPlotPosition = lastpos(1:2);

    % for TC plot window
    elseif strcmp(ch.Config.sattype, 'tcplot')

        % already remove field
        ne_gcfg.cc = rmfield(ne_gcfg.cc, varargin{3});

        % re-draw slice plot
        if ne_gcfg.fcfg.page < 3
            ne_setslicepos;
        end

        % return
        return;

    % for plot window
    elseif strcmp(ch.Config.sattype, 'mdmcondavg')

        % free MDM and VOI
        if isxff(ch.Config.mdm, true)
            ch.Config.mdm.ClearObject;
        end
        if isxff(ch.Config.voi, true)
            ch.Config.voi.ClearObject;
        end

        % record position
        ne_gcfg.c.ini.Children.MDMCondAvgPosition = lastpos(1:2);

    % for multi-image window
    elseif strcmp(ch.Config.sattype, 'multiimage')

        % delete Mixer window as well
        ch.Mixer.Delete;
        
        % then the image object
        if isxff(ch.Config.imobj, true)
            imobj = cc.Config.imobj;
            if ~isempty(imobj.FilenameOnDisk) && ...
                isfield(imobj.RunTimeVars, 'AutoSave') && ...
                islogical(imobj.RunTimeVars.AutoSave) && ...
                numel(imobj.RunTimeVars.AutoSave) == 1 && ...
                imobj.RunTimeVars.AutoSave
                cc.Config.imobj.SaveRunTimeVars;
            end
            cc.Config.imobj.ClearObject;
        end
        
        % and transio
        try
            delete(ch.Config.ti);
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
end
ne_gcfg.cc = rmfield(ne_gcfg.cc, varargin{3});
