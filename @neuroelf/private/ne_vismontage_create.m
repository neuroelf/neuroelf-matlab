% FUNCTION ne_vismontage_create: create montage image
function [m, ax, malp] = ne_vismontage_create(varargin)

% Version:  v1.0
% Build:    16011412
% Date:     Jan-14 2016, 12:42 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2016, Jochen Weber
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

% only allow one concurrent call
if any(strcmp('vismontagecreate', ne_gcfg.c.blockcb))
    return;
end
ne_gcfg.c.blockcb{end+1} = 'vismontagecreate';

% create options arguments
o = struct;

% get configuration
ipval = {'nearest', 'linear', 'cubic', 'lanczos3'};
cc = ne_gcfg.fcfg;
ch = ne_gcfg.h;
vm = cc.VisMontage;
tags = ch.VisMontage.h;

% anatomical transparency option and color
o.atrans = (tags.CB_vismontage_anatransp.Value > 0);
o.atranscol = lsqueeze(tags.BT_vismontage_anabackgc.CData(1, 1, :))';

% get blocking factor
blx = splittocellc( ...
    tags.DD_vismontage_layout.String{tags.DD_vismontage_layout.Value}, ' ');
o.blx = [str2double(blx{1}), str2double(blx{3})];

% and border size
o.brds = tags.DD_vismontage_imgborder.Value - 1;

% what column from the frame is direction and box?
switch (vm.dir)

    % sagittal slicing
    case {'sag'}

        % slice-through direction is X (1)
        o.drc = 1;

        % get stepping value
        o.drs = tags.DD_vismontage_xstep.Value;

    % for other directions
    case {'cor'}
        o.drc = 2;
        o.drs = tags.DD_vismontage_ystep.Value;
    case {'tra'}
        o.drc = 3;
        o.drs = tags.DD_vismontage_zstep.Value;
end

% for single slice
if o.drs > 10

    % set step to max box
    o.drs = 256;
end

% output filename
o.filename = ddeblank(tags.ED_vismontage_filename.String);

% flipping options (slice order, X-flipping)
o.flp = (tags.CB_vismontage_flipord.Value > 0);
o.flx = ((tags.CB_vismontage_flipX.Value > 0) && (o.drc > 1));

% font settings
o.fontcolor = lsqueeze(double(tags.BT_vismontage_fontcolor.CData(1, 1, :)))';
if any(o.fontcolor > 1)
    o.fontcolor = (1 / 255) .* o.fontcolor;
end
o.fontname = tags.DD_vismontage_fontname.String{tags.DD_vismontage_fontname.Value};
o.fontsize = str2double(strrep( ...
    tags.DD_vismontage_fontsize.String{tags.DD_vismontage_fontsize.Value}, 'pt', ''));

% get frame (with/without square/cube patch)
o.frame = vm.frame;

% copy of figure handle for visibility
o.hFig = vm.hFig;

% get interpolation technique
o.imeth = ipval{tags.DD_vismontage_interpst.Value};
o.imetha = ipval{tags.DD_vismontage_interpa.Value};

% joining stats maps
o.join = cc.join;

% get pixel per voxel
ppv = [0.5, 1, 2, 3, 4, 6, 8];
o.ppv = ppv(tags.DD_vismontage_pixpvox.Value);

% show-in-figure option
o.showinfig = (tags.RB_vismontage_showinfig.Value > 0);

% put on slice coordinates
o.slcoord = (tags.CB_vismontage_slcoord.Value > 0);

% get objects and index into maps
o.slvar = cc.SliceVar;
o.stvar = cc.StatsVar;
o.stvix = cc.StatsVarIdx;
if isempty(o.stvix)
    o.stalp = 1;
    o.stthr = [99, 999];
else
    stmap = o.stvar.Map(o.stvix(1));
    o.stalp = stmap.TransColorFactor;
    o.stthr = [stmap.LowerThreshold, stmap.UpperThreshold];
end

% show-while-slicing option
o.sws = (tags.CB_vismontage_dspslices.Value > 0);

% tempvol
o.tpvol = round(ch.Coord.TempSlider.Value);

% make the call
try
    % echo
    if ne_gcfg.c.echo
        oopts = o;
        if ~isxff(o.slvar, true) || ...
            isempty(o.slvar.FilenameOnDisk(true))
            oopts.slvar = '<VMR>';
        end
        if ~isxff(o.stvar, true) || ...
            isempty(o.stvar.FilenameOnDisk(true))
            oopts.stvar = '<VMP>';
        end
        ne_echo({'neuroelf_gui(''vismontage_create_ex'', %s);', ...
            any2ascii(oopts)});
    end
    [m, ax, malp] = ne_vismontage_create_ex(0, 0, o);
catch ne_eo;
    ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'vismontagecreate')) = [];
    rethrow(ne_eo);
end

% allow further montage creations
ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'vismontagecreate')) = [];
