% FUNCTION ne_srf_reload: reload currently selected SurfVar
function varargout = ne_srf_reload(varargin)

% Version:  v1.1
% Build:    16042117
% Date:     Apr-21 2016, 5:03 PM EST
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
cc = ne_gcfg.fcfg;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% get/check current SRF
srf = cc.SurfVar;
if numel(srf) ~= 1 || ...
   ~isxff(srf, 'srf')
    return;
end

% reload impossible
if isempty(srf.FilenameOnDisk)
    uiwait(warndlg('Surface not yet saved, cannot be reloaded.', ...
        'NeuroElf - info', 'modal'));
    return;
end

% echo
if ne_gcfg.c.echo
    ne_echo('srf', 'ReloadFromDisk');
end

% reload from disk
srf.ReloadFromDisk;

% and update colors
btc_meshcolor(srf, true);

% and update
srf.SetHandle('VertexCoordinateTal', []);
ne_setsurfpos(0, 0, 'upshape');
