% FUNCTION ne_vismontage_updateui: update UI based on callback
function ne_vismontage_updateui(varargin)

% Version:  v0.9b
% Build:    11082914
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

% get tags
vm = ne_gcfg.h.VisMontage;
fr = ne_gcfg.fcfg.VisMontage.frame;
tags = vm.h;

% depending on what to do
switch (lower(varargin{3}))

    % set to MNI brain box
    case {'brain'}

        % set pre-specified frame
        switch (tags.DD_vismontage_brainbox.Value)

            % full box
            case {1}
                ne_gcfg.fcfg.VisMontage.frame = ...
                    [128, 128, 128; -127.9999, -127.9999, -127.9999];

            % stats box
            case {2}
                try
                    nfr = 128 - ne_gcfg.fcfg.StatsVar.BoundingBox.BBox;
                    nfr = [nfr(1, [3, 1, 2]); nfr(2, [3, 1, 2]) + 0.01];
                    ne_gcfg.fcfg.VisMontage.frame = nfr;
                catch ne_eo;
                    ne_gcfg.c.lasterr = ne_eo;
                end

            % MNI brain
            case {3}
                ne_gcfg.fcfg.VisMontage.frame = ...
                    [78, 82, 80; -77.99, -109.99, -63.99];

            % AFNI brain
            case {4}
                ne_gcfg.fcfg.VisMontage.frame = ...
                    [72, 90, 64; -71.99, -89.99, -79.99];
        end
        ofr = ne_gcfg.fcfg.VisMontage.frame;

        % set text controls
        tags.ED_vismontage_xfrom.String = sprintf('%d', round(ofr(2, 1)));
        tags.ED_vismontage_xto.String = sprintf('%d', ofr(1, 1));
        tags.ED_vismontage_yfrom.String = sprintf('%d', round(ofr(2, 2)));
        tags.ED_vismontage_yto.String = sprintf('%d', ofr(1, 2));
        tags.ED_vismontage_zfrom.String = sprintf('%d', round(ofr(2, 3)));
        tags.ED_vismontage_zto.String = sprintf('%d', ofr(1, 3));

        % make sure the layout control is updated
        varargin{3} = 'dir';

    % change direction
    case {'dir'}

        % update the direction (used for OBJ::SliceToTransimg)
        ne_gcfg.fcfg.VisMontage.dir = ...
            tags.DD_vismontage_dir.String{tags.DD_vismontage_dir.Value}(1:3);

        % set the correct step control enabled
        tags.DD_vismontage_xstep.Enable = 'off';
        tags.DD_vismontage_ystep.Enable = 'off';
        tags.DD_vismontage_zstep.Enable = 'off';
        switch (tags.DD_vismontage_dir.Value)
            case {1}
                tags.DD_vismontage_xstep.Enable = 'on';
            case {2}
                tags.DD_vismontage_ystep.Enable = 'on';
            case {3}
                tags.DD_vismontage_zstep.Enable = 'on';
        end

    % change font color
    case {'fontcolor'}

        % old color
        ocol = lsqueeze(double(tags.BT_vismontage_fontcolor.CData(1, 1, :)))';
        if any(ocol > 1)
            ocol = (1 / 255) .* ocol;
        end
        ocol = min(255, max(0, round(255 .* ocol)));
        ncol = min(255, max(0, colorpicker(ocol, 'Slice label font color')));
        if any(ncol ~= ocol)
            tags.BT_vismontage_fontcolor.CData = ...
                repmat(reshape(uint8(ncol), [1, 1, 3]), 12, 12);
            drawnow;
        end

    % update from/to values
    case {'xfrom'}

        % get text value (with error allowance)
        try
            val = min(127, max(-128, round(str2double( ...
                tags.ED_vismontage_xfrom.String))));
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            val = -128;
        end

        % make sure it's not beyond other bound
        val = min(val, fr(1, 1));

        % update in frame
        ne_gcfg.fcfg.VisMontage.frame(2, 1) = val + 0.01;

        % and set again (rounding/error control)
        tags.ED_vismontage_xfrom.String = sprintf('%d', val);

    % for other controls
    case {'xto'}
        try
            val = min(128, max(-128, round(str2double( ...
                tags.ED_vismontage_xto.String))));
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            val = 128;
        end
        val = max(val, round(fr(2, 1)));
        ne_gcfg.fcfg.VisMontage.frame(1, 1) = val;
        tags.ED_vismontage_xto.String = sprintf('%d', val);
    case {'yfrom'}
        try
            val = min(127, max(-128, round(str2double( ...
                tags.ED_vismontage_yfrom.String))));
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            val = -128;
        end
        val = min(val, fr(1, 2));
        ne_gcfg.fcfg.VisMontage.frame(2, 2) = val + 0.01;
        tags.ED_vismontage_yfrom.String = sprintf('%d', val);
    case {'yto'}
        try
            val = min(128, max(-128, round(str2double( ...
                tags.ED_vismontage_yto.String))));
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            val = 128;
        end
        val = max(val, round(fr(2, 2)));
        ne_gcfg.fcfg.VisMontage.frame(1, 2) = val;
        tags.ED_vismontage_yto.String = sprintf('%d', val);
    case {'zfrom'}
        try
            val = min(127, max(-128, round(str2double( ...
                tags.ED_vismontage_zfrom.String))));
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            val = -128;
        end
        val = min(val, fr(1, 3));
        ne_gcfg.fcfg.VisMontage.frame(2, 3) = val + 0.01;
        tags.ED_vismontage_zfrom.String = sprintf('%d', val);
    case {'zto'}
        try
            val = min(128, max(-128, round(str2double( ...
                tags.ED_vismontage_zto.String))));
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            val = 128;
        end
        val = max(val, round(fr(2, 3)));
        ne_gcfg.fcfg.VisMontage.frame(1, 3) = val;
        tags.ED_vismontage_zto.String = sprintf('%d', val);

    % dis-/enable filename
    case {'show'}
        tags.ED_vismontage_filename.Enable = 'off';
        tags.CB_vismontage_slcoord.Enable = 'on';
        vm.VMFig.RadioGroupSetOne('VisMOut', 1);
    case {'write'}
        tags.CB_vismontage_slcoord.Value = 0;
        tags.CB_vismontage_slcoord.Enable = 'off';
        tags.ED_vismontage_filename.Enable = 'on';
        vm.VMFig.RadioGroupSetOne('VisMOut', 2);
end

% for frame update
if any('dxyz' == lower(varargin{3}(1)))

    % get new frame (rounded)
    fr = round(ne_gcfg.fcfg.VisMontage.frame);

    % depending on direction
    drc = tags.DD_vismontage_dir.Value;
    switch (drc)

        % X-dir
        case {1}

            % update layout if "dir" or an "x" -control
            upl = 'dx';

            % get stepsize
            drs = tags.DD_vismontage_xstep.Value;

        % for other directions
        case {2}
            upl = 'dy';
            drs = tags.DD_vismontage_ystep.Value;
        case {3}
            upl = 'dz';
            drs = tags.DD_vismontage_zstep.Value;
    end

    % for single slice
    if drs > 10

        % set step to max box
        drs = 256;
    end

    % if not the same as currently selection direction, return
    if ~any(upl == lower(varargin{3}(1)))
        return;
    end

    % get slicing positions
    slp = fr(2, drc):drs:fr(1, drc);
    slp(slp < -127 | slp > 128) = [];

    % compute possible layouts -> get number of slices
    nslp = numel(slp);

    % compute useful number of columns
    smls = unique([-2; -1; 0; 1; 2; 3; 4] + round(nslp .^ 0.4142));

    % and make sure each number is valid
    smls(smls < 1 | smls > nslp) = [];

    % create corresponding layout strings cell
    sstr = cell(numel(smls), 1);

    % and then compute the matching number of rows
    smls(:, 2) = ceil(nslp ./ smls);

    % and fill in strings
    for sc = 1:numel(sstr)
        sstr{sc} = sprintf('%d x %d', smls(sc, :));
    end

    % put into layout dropdown
    tags.DD_vismontage_layout.String = sstr;
    tags.DD_vismontage_layout.Value = ...
        min(numel(sstr), round(0.5 * numel(sstr)));
end
