% FUNCTION ne_fcprepro: preprocess VTCs for functional connectivity
function varargout = ne_fcprepro(varargin)

% Version:  v0.9c
% Build:    12051012
% Date:     May-10 2012, 12:42 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2012, Jochen Weber
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

% check slicevar
slvar = cc.SliceVar;
if numel(slvar) ~= 1 || ...
   ~isxff(slvar, {'fmr', 'hdr', 'vtc'}) || ...
    slvar.NrOfVolumes < 2
    return;
end
vols = slvar.NrOfVolumes;

% additional arguments given?
if nargin > 2 && ...
    ischar(varargin{3}) && ...
   ~isempty(varargin{3})
    outfile = varargin{3}(:)';
else
    [outfile, outpath] = uiputfile({'*.avi', 'AVI movie files (*.avi)'}, ...
        'Save new movie as...');
    if isequal(outfile, 0) || ...
        isequal(outpath, 0) || ...
        isempty(outfile)
        return;
    end
    outfile = [outpath '/' outfile];
end

% create AVI object
aviobj = avifile(outfile);

% undock current view
undocked = ne_undock;

% disallow resize
undocked.Resize = 'off';

% get handle and tag
uh = undocked.MLHandle;
ut = undocked.Tag(1:8);

% get a frame (for size estimate)
drawnow;
uf = getframe(uh);
ufd = uf.cdata;
ufs = size(ufd);
ufts = 16 .* ceil((1 / 16) .* ufs);
ufc = 1 + floor(0.5 .* (ufts - ufs));
ufd(1:ufts(1), 1:ufts(2), :) = 0;

% loop over volumes
for vc = 1:vols

    % set TempSlider
    ne_gcfg.cc.(ut).Coord.TempSlider.Value = vc;

    % update window and get frame
    ne_setsatslicepos(0, 0, ut);
    drawnow;
    uf = getframe(uh);

    % patch data (to multiple of 16)
    ufd(ufc(1):ufc(1)+ufs(1)-1, ufc(2):ufc(2)+ufs(2)-1, :) = uf.cdata;
    uf.cdata = ufd;

    % grab frame
    aviobj = addframe(aviobj, uf);
end

% save file
aviobj = close(aviobj);

% close window
ne_closesatwindow(0, 0, ut);

% return
if nargout > 0
    varargout{1} = aviobj;
end
