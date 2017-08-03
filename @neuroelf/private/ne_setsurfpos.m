function varargout = ne_setsurfpos(varargin)
% ne_setsurfpos  - update surface position and scenery
%
% FORMAT:       ne_setsurfpos([SRC, EVT, window, viewpt])
%
% Input fields:
%
%       SRC, EVT    Matlab handle callback inputs (discarded)
%       window      window specifier (used to target satellite windows)
%       viewpt      1xV cell array with order
%        {1}        1x1 double x-angle
%        {2}        1x1 double y-angle
%        {3}        1x3 translation (first value is discarded)
%        {4}        1x1 zoom
%        {5}        1x1 time (morph) index, starting at 0
%
% No output fields. (will be set to [])
%
% Example:
%
%     ne_setsurfpos(0, 0, 'BS123456', {180, 0, [0, -20, 0], 1.2, 1 / 3});
%
%     this sets the surface viewpoint in satellite window with ID
%     'BS123456' to the standard position (left-hemisphere) with a
%     zoom factor of 1.2 and a morph of 1/3 along the way between
%     the surfaces and the morphing targets

% Version:  v1.1
% Build:    16052714
% Date:     May-27 2016, 2:27 PM EST
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

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% get handles
if nargin < 3 || ~ischar(varargin{3}) || ~isfield(ne_gcfg.cc, varargin{3})
    fromroot = true;
    ch = ne_gcfg.h;
    pc = ne_gcfg.fcfg;

    % get mouse position
    nPos = get(ch.MainFigMLH, 'CurrentPoint');
    fPos = ne_gcfg.fcfg.surfpos;
    scn = mlhandle(ch.Scenery);
    sci = get(scn, 'Value');
    scu = get(scn, 'UserData');
else
    fromroot = false;
    sfid = varargin{3}(:)';
    ch = ne_gcfg.cc.(sfid);
    pc = ch.Config;
    nPos = ch.Satellite.CurrentPoint;
    fPos = ch.Config.surfpos;
    scn = ch.Scenery;
    sci = scn.Value;
    scu = scn.UserData;
end

% recoloring requested
recol = (~fromroot && nargin > 3 && isa(varargin{4}, 'double') && isequal(varargin{4}, 1));

% update of coordinates (shape) requested
upshape = (nargin > 2 && ischar(varargin{end}) && ~isempty(varargin{end}) && strcmpi(varargin{end}, 'upshape'));

% override settings
if nargin > 3 && iscell(varargin{4}) && ~isempty(varargin{4})

    % angle-x
    if isa(varargin{4}{1}, 'double') && numel(varargin{4}{1}) == 1 && ...
       ~isinf(varargin{4}{1}) && ~isnan(varargin{4}{1})

        % override
        if fromroot
            ne_gcfg.fcfg.srfcfg.anglex = mod(varargin{4}{1}, 360);
        else
            ne_gcfg.cc.(sfid).Config.srfcfg.anglex = mod(varargin{4}{1}, 360);
        end
    end

    % angle-y
    if numel(varargin{4}) > 1 && isa(varargin{4}{2}, 'double') && ...
        numel(varargin{4}{2}) == 1 && ~isinf(varargin{4}{2}) && ~isnan(varargin{4}{2})

        % override
        if fromroot
            ne_gcfg.fcfg.srfcfg.angley = min(90, max(-90, varargin{4}{2}));
        else
            ne_gcfg.cc.(sfid).Config.srfcfg.angley = min(90, max(-90, varargin{4}{2}));
        end
    end

    % translation
    if numel(varargin{4}) > 2 && isa(varargin{4}{3}, 'double') && ...
        any(numel(varargin{4}{3}) == [2, 3]) && ...
       ~any(isinf(varargin{4}{3}) | isnan(varargin{4}{3}))

        % override
        if fromroot
            ne_gcfg.fcfg.srfcfg.trans(1:numel(varargin{4}{3})) = ...
                min(256, max(-256, varargin{4}{3}(:)'));
        else
            ne_gcfg.cc.(sfid).Config.srfcfg.trans(1:numel(varargin{4}{3})) = ...
                min(256, max(-256, varargin{4}{3}(:)'));
        end
    end

    % zoom
    if numel(varargin{4}) > 3 && isa(varargin{4}{4}, 'double') && numel(varargin{4}{4}) == 1 && ...
       ~isinf(varargin{4}{4}) && ~isnan(varargin{4}{4}) && varargin{4}{4} > 0

        % override
        if fromroot
            ne_gcfg.fcfg.srfcfg.zoom = min(10, max(0.1, varargin{4}{4}));
        else
            ne_gcfg.cc.(sfid).Config.srfcfg.zoom = min(10, max(0.1, varargin{4}{4}));
        end
    end

    % time
    if numel(varargin{4}) > 4 && isa(varargin{4}{5}, 'double') && numel(varargin{4}{5}) == 1 && ...
       ~isinf(varargin{4}{5}) && ~isnan(varargin{4}{5}) && varargin{4}{5} >= 0
    
        % override
        if pc.srfcfg.time ~= varargin{4}{5}
            recol = true;
            upshape = true;
        end
        if fromroot
            ne_gcfg.fcfg.srfcfg.time = varargin{4}{5};
        else
            ne_gcfg.cc.(sfid).Config.srfcfg.time = varargin{4}{5};
        end
    end

    % re-read pc
    if fromroot
        pc = ne_gcfg.fcfg;
    else
        pc = ne_gcfg.cc.(sfid).Config;
    end
end

% position of controls on which clicks allows updating of position
cc = pc.srfcfg;

% do we need a hit-test
if ~upshape && any(ne_gcfg.c.btdown == gcbf) && isempty(pc.mpos.ddat)

    % make hit-test
    cobj = findfirst(fPos(:, 1) <= nPos(1) & fPos(:, 2) <= nPos(2) & ...
        fPos(:, 3) >  nPos(1) & fPos(:, 4) >  nPos(2));

    % object not hit?
    if isempty(cobj)

        % then pretend the button wasn't hit after all
        ne_gcfg.c.btdown = [];
        return;
    end

    % re-store config (as down-data)
    if fromroot
        ne_gcfg.fcfg.mpos.ddat = {3, nPos, cc};
    else
        ne_gcfg.cc.(sfid).Config.mpos.ddat = {3, nPos, cc};
    end

% we still need to update the position
elseif ~upshape && any(ne_gcfg.c.btdown == gcbf)

    % get original position and configuration
    oPos = pc.mpos.ddat{2};
    occ = pc.mpos.ddat{3};

    % depending on modifiers hit
    if isempty(pc.mpos.mods)

        % compute new angles
        cc.anglex = mod(occ.anglex + oPos(1) - nPos(1), 360);
        cc.angley = min(90, max(-90, occ.angley + 0.5 * (oPos(2) - nPos(2))));

    % shifting (translation)
    elseif numel(pc.mpos.mods) == 1 && strcmpi(pc.mpos.mods{1}, 'shift')

        % compute new translation
        cc.trans = min(256, max(-256, occ.trans + [0, nPos - oPos]));

    % zooming
    elseif numel(pc.mpos.mods) == 1 && strcmpi(pc.mpos.mods{1}, 'alt')

        % compute new translation
        cc.zoom = min(5, max(0.2, occ.zoom * (1.01 ^ (round(oPos(2) - nPos(2))))));

    % time point
    elseif numel(pc.mpos.mods) == 1 && strcmpi(pc.mpos.mods{1}, 'control')

        % update shape(s)
        upshape = true;

        % compute new time (for now, compute up to 300s, need config!)
        cc.time = round(min(1800, max(0, occ.time * 6 + nPos(1) - oPos(1)))) / 6;
        if pc.srfcfg.time ~= cc.time
            recol = true;
            upshape = true;
        end

        % for linked, take all
        if ne_gcfg.c.linked && (nargin < 4 || ~ischar(varargin{4}) || ...
            ~strcmpi(varargin{4}(:)', 'linkup'))
            ne_gcfg.fcfg.srfcfg.time = cc.time;
        elseif nargin > 3 && ischar(varargin{4}) && strcmpi(varargin{4}(:)', 'linkup')
            cc.time = ne_gcfg.fcfg.srfcfg.time;
        end

        % for each surface
        for sc = 1:numel(sci)

            % get handles
            scho = scu{sci(sc), 4};
            if ~isxff(scho, true)
                continue;
            end
            sch = handles(scho);
            if ~isfield(sch, 'Stats') || ~iscell(sch.Stats) || ...
                isempty(sch.Stats) || ~isxff(sch.Stats{1}, 'mtc')
                continue;
            end
            sch = sch.Stats{1};
            schrtv = sch.RunTimeVars;

            % compute new SubVolMap index
            if isfield(schrtv, 'AvgMTC') && schrtv.AvgMTC
                schmintime = schrtv.AvgWindowFrom;
                schstptime = schrtv.AvgWindowStep;
                sch.RunTimeVars.SubMapVol = ...
                    1 + (1000 * cc.time - schmintime) / schstptime;
            else
                sch.RunTimeVars.SubMapVol = 1 + 1000 * cc.time / sch.TR;
            end

            % update coloring
            btc_meshcolor(scho);
        end
    end

    % re-store config
    if fromroot
        ne_gcfg.fcfg.srfcfg = cc;
    else
        ne_gcfg.cc.(sfid).Config.srfcfg = cc;
    end
end

% linked
if ne_gcfg.c.linked && ...
   (nargin < 4 || ...
    ~ischar(varargin{4}) || ...
    ~strcmpi(varargin{4}(:)', 'linkup'))

    % temporarily set btdown to 0
    btd = ne_gcfg.c.btdown;
    ne_gcfg.c.btdown = 0;

    % update main window first
    try
        ne_gcfg.fcfg.srfcfg = cc;
        ne_setsurfpos(0, 0, [], 'linkup');

        % iterate over all satellite windows
        ccs = fieldnames(ne_gcfg.cc);
        ccs = ccs(~cellfun('isempty', regexp(ccs, '^BS')));
        for ccc = 1:numel(ccs)

            % if window is of type surf
            if strcmpi(ne_gcfg.cc.(ccs{ccc}).Config.sattype, 'surf')

                % update as well
                ne_gcfg.cc.(ccs{ccc}).Config.srfcfg = cc;
                ne_setsurfpos(0, 0, ccs{ccc}, 'linkup');
            end
        end
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end

    % return state
    ne_gcfg.c.btdown = btd;

    % return early!
    return;
end

% we updated
ne_gcfg.c.lastupd = now;

% viewpoint string
viewp = sprintf('azi: %d, zen: %.1f / xt: %d, yt: %d / zoom: %.3f (time: %.2fs)', ...
    cc.anglex, cc.angley, cc.trans(2), cc.trans(3), cc.zoom, cc.time);

% if not from root
if ~fromroot

    % iterate over visible surfaces
    if recol || upshape
        for uc = sci(:)'

            % and try to update vertices and normals with new config
            try
                singlesrf = scu{uc, 4};
                if recol
                    btc_meshcolor(singlesrf, true, sfid, scu{uc, 5}, true);
                end
                if upshape
                    [p, pn] = btc_meshcn(singlesrf, cc, ...
                        ~strcmpi(get(scu{uc, 5}, 'FaceColor'), 'none') || ...
                        ~strcmpi(pc.renderer, 'opengl'));
                    if ~isempty(singlesrf.TriangleVertex)
                        set(scu{uc, 5}, 'Vertices', p, 'VertexNormals', pn);
                    else
                        set(scu{uc, 5}, 'Vertices', p);
                    end
                end
            catch ne_eo;
                ne_gcfg.c.lasterr = ne_eo;
            end
        end
    end

    % general transform
    set(ch.SurfaceTransform, 'Matrix', btc_meshtrf(cc));

    % any plots on surface?
    pc = get(ch.Surface, 'Children');
    pc = pc(findfirst(strcmpi(get(pc, 'Type'), 'hggroup')));
    for pcc = 1:numel(pc)

        % check for valid UserData
        pcu = get(pc(pcc), 'UserData');
        if isstruct(pcu) && ...
            numel(pcu) == 1 && ...
            isfield(pcu, 'pcrd')

            % compute new coordinates
            pcrd = btc_coordc(pcu.pcrd, cc);
            if pcu.hover > 0
                pcrd(:, 1) = pcrd(:, 1) + pcu.hover;
            elseif pcu.hover < 0
                pcrd(:, 1) = -pcu.hover;
            end

            % set initial values
            set(pc(pcc), ...
                'CData', pcu.color, ...
                'XData', pcrd(:, 1), ...
                'YData', pcrd(:, 2), ...
                'ZData', pcrd(:, 3));
            if ~isempty(pcu.symsize)
                set(pc(pcc), 'SizeData', pcu.symsize);
            end
            if ~isempty(pcu.bcolor)
                set(pc(pcc), 'MarkerEdgeColor', pcu.bcolor);
            end

            % any texts
            if ~isnan(pcu.th(1))

                % shift positions
                if numel(pcu.labshift) == 1
                    labshift = pcrd(:, 2:3);
                    labshift = pcu.labshift .* (labshift ./ ...
                        repmat(sqrt(sum(labshift .* labshift, 2)), 1, 2));
                    pcrd(:, 2:3) = pcrd(:, 2:3) + labshift;
                else
                    pcrd(:, 2) = pcrd(:, 2) + pcu.labshift(1);
                    pcrd(:, 3) = pcrd(:, 3) + pcu.labshift(2);
                end
                for cpcc = 1:numel(pcu.th)
                    set(pcu.th(cpcc), ...
                        'Position', pcrd(cpcc, 1:3));
                end
            end
        end
    end

    % update other patches
    pc = get(ch.Surface, 'Children');
    pc = pc(strcmpi(get(pc, 'Type'), 'patch'));
    pc = setdiff(pc(:), cat(1, scu{sci, 5}, ch.SurfaceStatsBar));
    for pcc = 1:numel(pc)

        % check for valid UserData
        pcu = get(pc(pcc), 'UserData');
        if isstruct(pcu) && ...
            numel(pcu) == 1 && ...
            isfield(pcu, 'pcrd')

            % compute new coordinates
            pcrd = btc_coordc(pcu.pcrd, cc);
            set(pc(pcc), 'Vertices', pcrd(:, 1:3));
        end
    end

    % update name
    ch.Satellite.Name = sprintf('NeuroElf - satellite (%s) - %s', sfid, viewp);

    % return early
    return;
end

% update viewpoint
ch.SceneryViewPoint.String = viewp;

% make button available
if numel(sci) == 1 && ...
   ~isempty(scu)
    ch.SceneryProps.Enable = 'on';
else
    ch.SceneryProps.Enable = 'off';
end

% get visible flag
scb = false(size(scu, 1), 1);
scb(sci) = true;

% re-color flag
if nargin > 2 && ...
    islogical(varargin{3}) && ...
    numel(varargin{3}) == 1
    recol = varargin{3};
elseif nargin > 2 && ...
    isequal(varargin{3}, 1)
    recol = true;
end

% still invisible unless on page 3
if pc.page ~= 3
    visflag = 'off';
else
    visflag = 'on';
end

% update transform
set(ch.SurfaceTransform, 'Matrix', btc_meshtrf(cc));

% re-set visibility
for uc = 1:numel(scb)
    try
        f = scu{uc, 4};
        fh = handles(f);
        srfh = fh.Surface;
        srfp = fh.SurfProps;
        srft = fh.SurfaceTransform;

        % and if visible
        if scb(uc)

            % re-color
            if recol
                if fromroot
                    btc_meshcolor(f, recol);
                else
                    btc_meshcolor(f, true, sfid, scu{uc, 5}, true);
                end
            end

            % get current and requested patch mode
            facecmode = strcmpi(get(srfh, 'FaceColor'), 'interp');
            facermode = (srfp{5} == 'f');

            % get coordinates and normals and alpha
            if upshape
                [p, pn] = btc_meshcn(f, cc, facermode || ...
                    ~strcmpi(pc.renderer, 'opengl'));
            end

            % change in mode
            if facecmode ~= facermode

                % switch to face-mode
                if facermode
                    set(srfh, 'FaceColor', 'interp', 'EdgeColor', 'none', ...
                        'LineStyle', 'none');

                % switch to wireframe
                else
                    set(srfh, 'FaceColor', 'none', 'EdgeColor', 'interp', ...
                        'LineStyle', '-');
                end
            end

            % update patch
            if upshape
                if ~isempty(f.TriangleVertex)
                    set(srfh, 'Vertices', p, 'VertexNormals', pn, ...
                        'Visible', visflag);
                else
                    set(srfh, 'Vertices', p, ...
                        'Visible', visflag);
                end
            end
        else
            set(srfh, 'Visible', 'off');
        end
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
end

% update selection
if nargin > 2 && ...
    isa(varargin{3}, 'double') && ...
    isequal(varargin{3}, 1) && ...
    fromroot

    % update bar status
    ne_setcsrfstatbars;

    % if only one file is shown
    if numel(sci) == 1

        % get the requested object
        singlesrf = scu{sci, 4};
        
        % set name
        [srfpath, ne_gcfg.c.title{2, 1}] = fileparts(singlesrf.FilenameOnDisk);
        if isempty(ne_gcfg.c.title{2, 1})
            ne_gcfg.c.title{2, 1} = 'unsaved';
        end
        
        % stats?
        srfhandles = singlesrf.Handles;
        if isfield(srfhandles, 'Stats') && ...
            iscell(srfhandles.Stats) && ...
            numel(srfhandles.Stats) == 2 && ...
            numel(srfhandles.Stats{1}) == 1 && ...
            isxff(srfhandles.Stats{1}, true)
            [smppath, ne_gcfg.c.title{2, 2}] = fileparts(srfhandles.Stats{1}.FilenameOnDisk);
            if isempty(ne_gcfg.c.title{2, 2})
                ne_gcfg.c.title{2, 2} = 'unsaved';
            end
            if numel(srfhandles.Stats{2}) == 1 && ...
                srfhandles.Stats{2} <= srfhandles.Stats{1}.NrOfMaps
                ne_gcfg.c.title{2, 3} = srfhandles.Stats{1}.Map(srfhandles.Stats{2}).Name;
            else
                ne_gcfg.c.title{2, 3} = sprintf('(%d maps)', numel(srfhandles.Stats{2}));
            end
        end

        % then get the userdata of the dropdown
        srfs = ch.SurfVar.UserData;
        for uc = 1:size(srfs, 1)
            if singlesrf == srfs{uc, 4}
                ch.SurfVar.Value = uc;
                ne_setcsrf(0, 0, false);
                break;
            end
        end
    else
        ne_gcfg.c.title(2, :) = {sprintf('%d surfaces', numel(sci)), '', ''};
    end
    
    % update name
    ne_updatename;
end

% remote
if ne_gcfg.c.remote

    % grab screenshot and write to images folder
    simg = getframe(ch.Surface);
    ifmt = lower(ne_gcfg.c.ini.Remote.ImageFormat);
    if strcmp(ifmt, 'jpg')
        iqual = {'Quality', ne_gcfg.c.ini.Remote.ImageJPGQuality};
    else
        iqual = {};
    end
    ipath = [neuroelf_path('remote') '/images'];
    try
        imwrite(simg.cdata, sprintf('%s/surface.%s', ipath, ifmt), iqual{:});
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
end
