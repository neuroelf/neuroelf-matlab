% FUNCTION ne_setcvar: set current SliceVar (get from control)
function varargout = ne_setcvar(varargin)

% Version:  v1.1
% Build:    16052611
% Date:     May-26 2016, 11:37 AM EST
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
ch = ne_gcfg.h;
cf = ch.MainFig;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% set initial visible groups
cf.SetGroupVisible('FMRMenu', 'off');
cf.SetGroupVisible('VMRMenu', 'off');
cf.SetGroupVisible('VTCMenu', 'off');
cf.SetGroupVisible('HDRMenu', 'off');
cf.SetGroupVisible('DrawMenu','off');
cf.SetGroupEnabled('SVIsVMR', 'off');
cf.SetGroupEnabled('VHasBVs', 'off');
cf.SetGroupEnabled('VIsTemp', 'off');
set(ch.TCPlotDiscards, 'Visible', 'off');

% check DTI dataset
if numel(ne_gcfg.fcfg.dti) ~= 1 || ...
   ~isxff(ne_gcfg.fcfg.dti, 'hdr')
    ch.DTIPlotFibers.Enable = 'off';
    ch.DTITrackFibers.Enable = 'off';
end

% and force some menu and button items
ch.Menu.LimitVMR.Enable = 'off';
ch.VMRShowV16.Enable = 'off';
ch.VMRShowV16.Value = 0;
showv16 = false;

% unset any SliceUnder var
ne_gcfg.fcfg.SliceUnder = ...
    struct('Filetype', 'NONE', 'RunTimeVars', struct('Trf', eye(4)));

% get data from control
hu = ch.SliceVar.UserData;
if nargin < 3 || ...
   ~isa(varargin{3}, 'double') || ...
    numel(varargin{3}) ~= 1 || ...
    isinf(varargin{3}) || ...
    isnan(varargin{3}) || ...
    varargin{3} < 1 || ...
    varargin{3} > size(hu, 1) || ...
    varargin{3} ~= fix(varargin{3})
    hi = ch.SliceVar.Value;
else
    hi = varargin{3};
    ch.SliceVar.Value = hi;
end

% preset title
ne_gcfg.c.title{1, 1} = ch.SliceVar.String{ch.SliceVar.Value};

% check validity
if hi < 1 || ...
    hi > size(hu, 1)
    ne_gcfg.fcfg.SliceVar = struct('Filetype', 'NONE');
    if isstruct(ne_gcfg.fcfg.Render) && ...
        isfield(ne_gcfg.fcfg.Render, 'hFig') && ...
        isxfigure(ne_gcfg.fcfg.Render.hFig, true)
        ne_render_closeui;
    end
    ne_updatename;
    return;
end

% vars available
if ~isempty(hu)

    % get xff object from UserData
    slvar = hu{hi, 4};
    slrtv = slvar.RunTimeVars;

    % and set certain groups
    cf.SetGroupEnabled('VarIsLoaded', 'on');

% no vars
else

    % default
    slvar = struct('Filetype', 'NONE');
    slrtv = struct;

    % and set certain groups
    cf.SetGroupEnabled('VarIsLoaded', 'off');
end
slvartyp = lower(slvar.Filetype);

% ensure RTV has some minimal fields
if ~isfield(slrtv, 'TrfPlus')
    slvar.RunTimeVars.TrfPlus = eye(4);
end
if ~isfield(slrtv, 'AlphaVolume')
    slvar.RunTimeVars.AlphaVolume = '';
end

% default cursor step
ne_gcfg.fcfg.cstep = 1;

% no painting
ne_gcfg.fcfg.paint.mode = 1;

% put into figure configuration (so control hasn't to be accessed again)
ne_gcfg.fcfg.SliceVar = slvar;

% get BV-style position (for controls and value sampling) -> 3D BV-types
if ~any(strcmp(slvartyp, {'dmr', 'fmr', 'hdr', 'head', 'nlf'}))

    % use bvcoordconv
    if ~isfield(slrtv, 'Trf')
        slvar.RunTimeVars.Trf = ...
            bvcoordconv(zeros(0, 3), 'tal2bvc', slvar.BoundingBox)';
    end

% otherwise
else

    % buffer Trf
    if ~isfield(slrtv, 'Trf')
        slvar.RunTimeVars.Trf = inv(slvar.CoordinateFrame(1).Trf)';
    end
end

% reset temp-slider (anyway)
ch.Coord.TempSlider.Value = 1;
ch.Coord.Temp.String = '1';

% for functional types
if any(strcmp(slvartyp, {'dmr', 'fmr', 'vtc'})) || ...
   (strcmp(slvartyp, 'hdr') && size(slvar.VoxelData, 4) > 1) || ...
   (strcmp(slvartyp, 'head') && numel(slvar.Brick) > 1)

    % enable the TC-plot
    ne_gcfg.fcfg.tcplot = true;

    % update for VTCs
    slvarnvol = slvar.NrOfVolumes;
    if strcmp(slvartyp, 'vtc')

        % also set VTC-discard image
        tcdisc = repmat(reshape(uint8([255, 0, 0]), [1, 1, 3]), 1, slvarnvol);
        tcdisv = zeros(1, slvarnvol);
        if isfield(slrtv, 'Discard') && ~isempty(slrtv.Discard)
            tcdisv(slrtv.Discard) = 0.5;
        end
        set(ch.TCPlotDiscards, 'CData', tcdisc, 'AlphaData', tcdisv, ...
            'XData', [1, slvarnvol], 'YData', [-32768, 32768]);

        % instant correlations?
        if ne_gcfg.fcfg.instscorr
            ne_gcfg.fcfg.instscvar = slvar;
        end
    end

    % and the controls of it
    cf.SetGroupEnabled('VIsTemp', 'on');

    % and set a valid slider range
    ch.Coord.TempSlider.Max = max(1, slvarnvol);
    if slvar.NrOfVolumes > 1
        ch.Coord.TempSlider.SliderStep = ...
            max(eps, (1 / (slvar.NrOfVolumes - 1)) * [1, 10]);
    else
        ch.Coord.TempSlider.SliderStep = [1, 1];
    end

    % for DTI-based HDRs
    if strcmp(slvartyp, 'hdr') && ...
        isfield(slvar.RunTimeVars, 'BVals') && ...
        isfield(slvar.RunTimeVars, 'BVecs')
        cf.SetGroupEnabled('VHasBVs', 'on');
        ne_gcfg.fcfg.dti = slvar;
        ch.DTIPlotFibers.Enable = 'on';
        ch.DTITrackFibers.Enable = 'on';
    end

% otherwise
else

    % disable TC-plot
    ne_gcfg.fcfg.tcplot = false;

    % already disable display (before...)
    set(ch.TCPlot.Children, 'Visible', 'off');
    ch.TCPlot.Visible = 'off';
    cf.SetGroupEnabled('TCPlot', 'off');
    ccf = fieldnames(ne_gcfg.cc);
    if ~isempty(ccf)
        ccf(cellfun('isempty', regexp(ccf, '^TC'))) = [];
    end
    if ~isempty(ccf) && numel(ne_gcfg.cc.(ccf{1})) == 1 && isstruct(ne_gcfg.cc.(ccf{1})) && ...
        isfield(ne_gcfg.cc.(ccf{1}), 'Satellite') && isxfigure(ne_gcfg.cc.(ccf{1}).Satellite, true)
        tcvar = ne_gcfg.cc.(ccf{1}).Config.tcvar;
        if ~isxff(tcvar, 'vtc') || ~isfield(tcvar.RunTimeVars, 'AvgVTC') || ...
           ~islogical(tcvar.RunTimeVars.AvgVTC) || ~tcvar.RunTimeVars.AvgVTC
            ne_gcfg.cc.(ccf{1}).Satellite.Visible = 'off';
        end
    end

    % for VMRs with 1x1x1 resolutin and Analyze: enable plotting commands
    if strcmp(slvartyp, 'vmr') || ...
        strcmp(slvartyp, 'hdr')
        cf.SetGroupVisible('DrawMenu', 'on');
        cf.SetGroupEnabled('SVIsVMR', 'on');

        % lower cursor step?
        if strcmp(slvartyp, 'vmr') && ...
            slvar.VoxResX == slvar.VoxResY && ...
            slvar.VoxResX == slvar.VoxResZ && ...
            slvar.VoxResX < 1
            ne_gcfg.fcfg.cstep = slvar.VoxResX;
        elseif strcmp(slvartyp, 'hdr') && ...
            max(slvar.CoordinateFrame.Resolution) < 1
            ne_gcfg.fcfg.cstep = 2 ^ round(log(mean(slvar.CoordinateFrame.Resolution)));
        end

        % with V16 data
        if strcmp(slvartyp, 'vmr')
            if ~isfield(slvar.RunTimeVars, 'ShowV16')
                slvar.RunTimeVars.ShowV16 = false;
            end
            if isequal(size(slvar.VMRData), size(slvar.VMRData16))

                % enable LimitVMR menu item and V16 button
                ch.Menu.LimitVMR.Enable = 'on';
                ch.VMRShowV16.Enable = 'on';

                % check value
                showv16 = slvar.RunTimeVars.ShowV16;
                ne_gcfg.fcfg.showv16 = showv16;
                if showv16
                    ch.VMRShowV16.Value = 1;
                end
            end
        end
    else
        cf.SetGroupEnabled('SVIsVMR', 'off');
    end
end

% set correct group visible
if any(strcmp(slvartyp, {'fmr', 'vmr', 'vtc', 'hdr'}))
    cf.SetGroupVisible([upper(slvartyp) 'Menu'], 'on');
end

% update histogram of scaling window
if ~isfield(slvar.RunTimeVars, 'ScalingHist') || ...
   (strcmp(slvartyp, 'vmr') && ...
    ~isempty(slvar.VMRData16) && ...
    ~isfield(slvar.RunTimeVars, 'ScalingHist16'))
    slvar.SetScalingWindow;
end
rtv = slvar.RunTimeVars;
if showv16
    swl = rtv.ScalingWindowLim16;
    swn = rtv.ScalingWindow16 - swl(1);
    shst = rtv.ScalingHist16(:);
else
    swl = rtv.ScalingWindowLim;
    swn = rtv.ScalingWindow - swl(1);
    shst = rtv.ScalingHist(:);
end
swld = abs(swl(2) - swl(1));
swn(1) = max(0, swn(1));
swnd = abs(swn(2) - swn(1));
shst = flexinterpn(shst, [Inf; 1; 1; numel(shst)], [0; 0.5; 1; 1; 1; 0.5; 0], 1);
sshst = sort(shst(shst > 0));
if numel(sshst) > 16
    shst = (0.5 / mean(sshst(ceil(0.5 * numel(sshst)):end - 8))) .* shst;
elseif ~isempty(sshst)
    shst = (1 / max(sshst)) .* shst;
else
    shst = zeros(256, 1);
end
set(ch.HistPlot, 'XData', 0.1 + min(0.8, 0.8 .* shst));
hval = swn ./ swld;
ne_gcfg.fcfg.histval = hval;
set(ch.HistLine1, 'YData', [hval(1); hval(1)]);
set(ch.HistLine2, 'YData', [hval(2); hval(2)]);
hid = (0:255)';
hid = uint8(min(255, max(0, round( ...
    (hid - (256 / swld) * swn(1)) * (swld / swnd)))));
set(ch.HistImage, 'CData', repmat(hid, [1, 16, 3]));

% underlay data?
slhnd = handles(slvar);
if isfield(slhnd, 'Underlay') && ...
    numel(slhnd.Underlay) == 1 && ...
    isxff(slhnd.Underlay, {'fmr', 'hdr', 'head', 'mgh', 'vmr', 'vtc'})
    ne_gcfg.fcfg.SliceUnder = slhnd.Underlay;
end

% show different page ?
if ne_gcfg.fcfg.page > 2 && ...
    ne_gcfg.fcfg.page ~= 4
    ne_showpage(0, 0, 1);
end

% update render page
if isstruct(ne_gcfg.fcfg.Render) && ...
    isfield(ne_gcfg.fcfg.Render, 'hFig') && ...
    isxfigure(ne_gcfg.fcfg.Render.hFig, true)

    % close UI?
    if ~isxff(slvar, {'hdr', 'head', 'vmr', 'vtc'})
        ne_render_closeui;
        ne_setslicepos;
        ne_updatename;
        return;
    end

    ne_gcfg.fcfg.Render.slvar = slvar;
    slavar = rtv.AlphaVolume;
    if ischar(slavar) && ...
        numel(slavar) == 24 && ...
        all((slavar(:) >= '0' & slavar(:) <= '9') | (lower(slavar(:)) >= 'a' & lower(slavar(:)) <= 'f'))
        try
            slavar = xff(lower(slavar(:)'));
            if ~isxff(slavar, {'hdr', 'head', 'vmr', 'vtc'})
                slavar = [];
            end
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            slavar = [];
        end
    else
        slavar = [];
    end
    slvarfile = slvar.FilenameOnDisk(2);
    if isempty(slvarfile)
        slvarfile = sprintf('<untitled.%d>', lower(slvar.Filetype));
    elseif numel(slvarfile) > 51
        slvarfile = [slvarfile(1:24) '...' slvarfile(end-23:end)];
    end
    slvars = ch.SliceVar.UserData;
    slvarsn = {'<same>'};
    slvarso = {[]};
    slvci = 1;
    for slvc = 1:size(slvars, 1)
        if isxff(slvars{slvc, 4}, true) && ...
            slvars{slvc, 4} ~= slvar
            slavarc = slvars{slvc, 4};
            slvarsn{end+1} = slvars{slvc, 1};
            slvarso{end+1} = slavarc;
            if slvars{slvc, 4} == slavar
                slvci = numel(slvarsn);
            end
            slavarh = handles(slavarc);
            if isfield(slavarh, 'RenderAVol')
                slavarc.DeleteHandle('RenderAVol');
            end
            if isfield(slavarh, 'RenderSVol')
                slavarc.DeleteHandle('RenderSVol');
            end
            if isfield(slavarh, 'RenderPreview')
                slavarc.DeleteHandle('RenderPreview');
            end
            if isfield(slavarh, 'RenderPView')
                slavarc.DeleteHandle('RenderPView');
            end
        end
    end
    hTag = ne_gcfg.fcfg.Render.hTag;
    hTag.ED_render_slvarfile.String = ['  ' slvarfile];
    hTag.DD_render_slavarfile.Value = 1;
    hTag.DD_render_slavarfile.String = slvarsn;
    hTag.DD_render_slavarfile.UserData = slvarso;
    if numel(slvarsn) > 1
        ne_gcfg.fcfg.Render.hFig.SetGroupEnabled('SlAlpha', 'on');
        hTag.DD_render_slavarfile.Value = slvci;
    else
        ne_gcfg.fcfg.Render.hFig.SetGroupEnabled('SlAlpha', 'off');
    end
end

% update the window
ne_setslicepos;
ne_updatename;

% also update render view
if ne_gcfg.fcfg.page == 4 && ...
    isstruct(ne_gcfg.fcfg.Render)
    ne_render_setview;
end
