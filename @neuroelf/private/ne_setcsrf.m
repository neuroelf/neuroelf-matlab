% PUBLIC FUNCTION ne_setcsrf: set current SurfVar (get from control)
function varargout = ne_setcsrf(varargin)

% Version:  v1.1
% Build:    16031521
% Date:     Mar-15 2016, 9:51 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2014, 2016, Jochen Weber
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
cf = ch.MainFig;

% default: hide SRFMenu
cf.SetGroupVisible('SRFMenu', 'off');

% valid input given
hi = [];
hu = ch.SurfVar.UserData;
ssp = true;
if nargin > 2

    % index
    if isa(varargin{3}, 'double') && ...
        numel(varargin{3}) == 1 && ...
        any(1:size(hu, 1) == varargin{3})
        hi = varargin{3};

    % xff
    elseif numel(varargin{3}) == 1 && ...
        isxff(varargin{3}, {'fsbf', 'srf'})

        % in list
        for ic = 1:size(hu, 1)
            if hu{ic, 4} == varargin{3}
                hi = ic;
                break;
            end
        end
        
    % setsurfpos
    elseif numel(varargin{3}) == 1 && ...
        islogical(varargin{3})
        ssp = false;
    end
end

% get data from control
if isempty(hi)
    hi = ch.SurfVar.Value;
end

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% check validity
if hi < 1 || ...
    hi > size(hu, 1)
    ne_gcfg.fcfg.SurfVar = struct('Filetype', 'NONE');
    return;
end

% vars available
if ~isempty(hu)

    % get xff object from UserData
    srfvar = hu{hi, 4};

    % and set certain groups
    cf.SetGroupEnabled('VarIsLoaded', 'on');

% no vars
else

    % default
    srfvar = struct('Filetype', 'NONE');

    % and set certain groups
    cf.SetGroupEnabled('VarIsLoaded', 'off');
end

% put into figure configuration (so control hasn't to be accessed again)
ne_gcfg.fcfg.SurfVar = srfvar;
if nargout > 0
    varargout{1} = srfvar;
end

% for real SRFs
if isxff(srfvar, true)

    % set correct group visible
    cf.SetGroupVisible('SRFMenu', 'on');

    % show different page ?
    if ne_gcfg.fcfg.page < 3
        ne_showpage(0, 0, 3);
    end

    % stats object linked?
    srfh = handles(srfvar);
    luo = srfh.Stats{1};
    lud = ch.SurfStatsVar.UserData;
    ssmi = false;
    if ~isempty(luo) && ...
        isxff(luo, true)

        % then lookup
        for slc = 1:size(lud, 1)
            if lud{slc, 4} == luo
                ch.SurfStatsVar.Value = slc;
                ssmi = true;
                break;
            end
        end

    % no object linked
    else

        % look through SurfStatsVars
        nrvert = srfvar.NrOfVertices;
        for slc = 1:size(lud, 1)
            if isxff(lud{slc, 4}, true) && ...
                lud{slc, 4}.NrOfVertices == nrvert && ...
               ~isempty(lud{slc, 4}.Map)
                ch.SurfStatsVar.Value = slc;
                srfvar.SetHandle('Stats', {lud{slc, 4}, 1});
                srfh.Stats = {lud{slc, 4}, 1};
                ssmi = true;
                break;
            end
        end
    end

    % update stats
    ne_setcsrfstats;
    if ssmi
        ch.SurfStatsVarMaps.Value = srfh.Stats{2};
        ne_setcsrfstatmap;
    end
end

% update the window
if ssp
    ne_setsurfpos;
end
