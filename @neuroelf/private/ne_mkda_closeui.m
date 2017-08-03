% FUNCTION ne_mkda_closeui: close MKDA window
function ne_mkda_closeui(varargin)

% Version:  v0.9c
% Build:    11120213
% Date:     Nov-10 2011, 10:32 AM EST
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

% disallow while being run
if any(strcmp(ne_gcfg.c.blockcb, 'mkda_run'))
    return;
end

% close all loaded PLPs
plpud = ch.PLPs.UserData;
plpst = ch.PLPs.String;
if ~iscell(plpst)
    plpst = cellstr(plpst);
end
for pc = size(plpud, 1):-1:1
    if numel(plpud{pc, 3}) == 1 && ...
        isxff(plpud{pc, 3}, 'plp')
        plp = plpud{pc, 3};
        if ~plp.RunTimeVars.Saved
            plpid = plp.FilenameOnDisk;
            if isempty(plpid)
                if isfield(plp.RunTimeVars, 'SourceFile') && ...
                    ischar(plp.RunTimeVars.SourceFile) && ...
                   ~isempty(plp.RunTimeVars.SourceFile)
                    plpid = plp.RunTimeVars.SourceFile;
                else
                    plpid = ['xffID #', plp.RunTimeVars.xffID];
                end
            end
            vans = questdlg(sprintf('PLP object %s unsaved. Save now?', plpid), ...
                'NeuroElf - user input', 'Yes', 'No', 'Cancel', 'Yes');
            if isempty(vans) || ...
               ~ischar(vans) || ...
                strcmpi(vans(:)', 'cancel')
                return;
            end
            if strcmpi(vans(:)', 'yes')
                plp.SaveAs;
            end
        end
        plp.ClearObject;
        plpud(pc, :) = [];
        plpst(pc) = [];
        ch.PLPs.UserData = plpud;
        if isempty(plpst)
            plpst = {'<no PLP loaded>'};
        end
        ch.PLPs.Value = numel(plpst);
        ch.PLPs.String = plpst;
        ne_mkda_setplp;
    end
end

% update last known position
ne_gcfg.c.ini.Children.MKDAPosition = hFig.Position(1:2);

% delete figure and remove from global struct
hFig.Delete;
ne_gcfg.fcfg.MKDA = [];
ne_gcfg.h.MKDA = [];

% release contrast manager call
ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'mkda_open')) = [];
