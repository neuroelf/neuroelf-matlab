% FUNCTION ne_vmp_createmsk: create MSK from current VMP
function varargout = ne_vmp_createmsk(varargin)

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

% only valid if single VMP
stvar = ne_gcfg.fcfg.StatsVar;
if numel(stvar) ~= 1 || ...
   ~isxff(stvar, 'vmp')
    return;
end

% get maps to create MSK from
mapsel = ch.StatsVarMaps.Value;
if isempty(mapsel)
    return;
end

% options
if numel(mapsel) > 1
    opts = inputdlg( ...
        {'How to combine maps: (and/mean/or)', ...
         'Mean threshold: (]0 .. 1])'}, 'NeuroElf - user input', 1, ...
        {'  and', '  0.5'});
    if ~iscell(opts) || ...
        numel(opts) ~= 2 || ...
       ~any(strcmpi(ddeblank(opts{1}), {'and', 'mean', 'or'}))
        return;
    end
    opts = struct('combtype', ddeblank(opts{1}), 'mthresh', str2double(opts{2}));
    if numel(opts.mthresh) ~= 1 || ...
        isinf(opts.mthresh) || ...
        isnan(opts.mthresh) || ...
        opts.mthresh <= 0 || ...
        opts.mthresh > 1
        return;
    end
else
    opts = struct;
end

% generate MSK
msk = cell(1, 1);
try
    msk{1} = stvar.CreateMSK(mapsel, opts);

    % echo
    if ne_gcfg.c.echo
        ne_echo('vmp', 'CreateMSK', mapsel, opts);
    end

    % call SaveAs
    msk{1}.SaveAs;
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
end
clearxffobjects(msk);
