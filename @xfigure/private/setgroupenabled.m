function varargout = setgroupenabled(xo, iStr, emode, rflag)
%XFIGURE::SETGROUPENABLED  Sets a group of elements (EGroup) enabled.
%   SETGROUPENABLED(FIG, GROUP) sets the group GROUP in figure FIG
%   enabled.
%
%   SETGROUPENABLED(FIG, GROUP, EMODE) sets the enabled state of the group
%   GROUP in figure FIG to EMODE (which can also be 1x1 numeric).

% Version:  v1.1
% Build:    16041813
% Date:     Apr-18 2016, 1:30 PM EST
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

% global function library
global ne_methods;

% set output
varargout = cell(1, nargout);

% only valid for figures
if numel(xo) ~= 1 || xo.T ~= 1
    error('neuroelf:xfigure:invalidObjectType', 'SetGroupEnabled is only valid for figures.');
end

% check group spec
if ~ischar(iStr) || isempty(iStr(:))
    error('neuroelf:xfigure:invalidGroupName', 'SetGroupEnabled requires a valid group name.');
end
dogroup = deblank(iStr(:)');
enamode = 'on';

% determine target mode
if nargin > 2 && ~isempty(emode)
    if (((isnumeric(emode) || islogical(emode)) && ~emode(1)) || ...
         (ischar(emode) && strcmpi(emode, 'off')))
        enamode = 'off';
    elseif isnumeric(emode) && (ishandle(emode(1)) || ishandle(-emode(1)))
        if emode < 0
            isnegated = true;
            emode = -emode(1);
        else
            isnegated = false;
            emode =  emode(1);
        end
        try
            tHnd = xfigure(emode);
            if isxfigure(tHnd, 1)
                vls = ne_methods.subget(tHnd.mhnd, {'Value', 'Min', 'Max'});
                if vls.Value ~= vls.Max && isnegated
                    enamode = 'off';
                end
            end
        catch xfigerror
            neuroelf_lasterr(xfigerror);
            return;
        end
    end
end

% all groups ?
if strcmpi(dogroup, 'all_groups')
    gnames = fieldnames(xo.X.egroups);
    for gc = 1:length(gnames)
        setgroupenabled(xo, gnames{gc}, enamode);
    end
    return;
elseif any(dogroup == ';' | dogroup == ',' | dogroup == ' ')
    gnames = ne_methods.splittocellc(dogroup, ';, ', true, true);
    for gc = 1:length(gnames)
        if isvarname(gnames{gc})
            setgroupenabled(xo, gnames{gc}, enamode);
        end
    end
    return;
elseif ~isvarname(deblank(dogroup))
    error('neuroelf:xfigure:invalidGroupName', 'SetGroupEnabled requires a valid group name.');
end

% single group
if isfield(xo.X.figprops.egroups, dogroup)
    hdls = xo.X.figprops.egroups.(dogroup);
    for hc = 1:numel(hdls)
        setprop(hdls(hc), 'Enable', enamode);
    end
end

% refresh graphics
if ((nargin > 2 && ischar(emode) && strcmpi(emode, 'refresh')) || ...
    (nargin > 3 && ischar(rflag) && strcmpi(rflag, 'refresh'))) 
    redrawfig(xo.H);
end
