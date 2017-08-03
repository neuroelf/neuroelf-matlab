function p = xprogress(p, pval, ptxt, pstate, pmin, pmax)
% xprogress  - show a progress bar
%
% FORMAT:       p = xprogress([opts]);
%         or    [p] = xprogress(p, pval, ptxt, pstate, pmin, pmax)
%         or    xprogress(p, cmd, cmdarg);
%
% Input fields:
%
%       opts        1x1 struct with optional settings
%        .colors    2x3 RGB for bar colors
%
%       p           xprogress object
%       pval        current progress value
%       ptxt        label (e.g. current progress step)
%       pstate      one of 'close', 'hidden', {'visible'}
%       pmin        minimum value, default: 0
%       pmax        maximum value, default: 1
%
%       cmd         either of
%                   'close' (use an empty 3rd arg)
%                   'setcaption'
%                   'setposition' (1x2 screen pixels)
%                   'setsize' (1x2 screen size)
%                   'settitle'
%                   'setvisible'
%       cmdarg      a valid argument for the command
%
% Note: to start using the class, create an empty progress object
%       which will have the state hidden
%
%       p = xprogress;
%
%       then call the class method to change the object
%
%       [p] = xprogress(p, ...);
%       xprogress(p, 'settitle', 'Preparing GLM...');

% Version:  v1.1
% Build:    16032821
% Date:     Mar-28 2016, 9:00 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
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

% persistent configuration
persistent xpconfig;
if isempty(xpconfig)
    xpconfig = struct('lastid', -99, 'lastpos', [-1, -1, -1, -1], 'lasttxt', '');
    nelf = neuroelf;
end

% requires neuroelf
global ne_methods;

% default constructor
if nargin < 1 || numel(p) ~= 1 || ~strcmpi(class(p), 'xprogress')

    % but only called for empty argument list
    if nargin > 0 && (~isstruct(p) || numel(p) ~= 1)
        error('neuroelf:xprogress:badCall', 'Bad call to xprogress.');
    end
    if nargin > 0
        opts = p;
    else
        opts = struct;
    end

    % MATLAB version
    mlversion = str2double(ne_methods.splittocell(regexprep(version, ' .*$', ''), '.'));
    
    % get ROOT object properties
    rp = get(0);
    if ~any(rp.ScreenSize > 1)
        error('neuroelf:xprogress:noScreen', 'No GUI object available.');
    end

    % generate figure object
    hFig = figure(...
        'CloseRequestFcn',  'try, set(gcbf, ''Visible'', ''off''); end', ...
        'HandleVisibility', 'off', ...
        'IntegerHandle',    'off', ...
        'Menubar',          'none', ...
        'Name',             'Progress...', ...
        'NumberTitle',      'off', ...
        'Toolbar',          'none', ...
        'Units',            'pixels', ...
        'Visible',          'off');
    set(hFig, 'CloseRequestFcn', 'try, set(gcbf, ''Visible'', ''off''); end');
    set(hFig, 'Visible',   'off');
    if all(xpconfig.lastpos == -1)
        set(hFig, 'Position', [rp.ScreenSize(3) / 2 - 240, rp.ScreenSize(4) / 2, 480, 36]);
    else
        set(hFig, 'Position', xpconfig.lastpos);
    end

    % configure progress bar objects
    numlines = 128;
	numcols  = fix((numlines + 1) / 2);
    iCap = 'Progress:';
    if isfield(opts, 'color') && ...
        isa(opts.color, 'double') && ...
        numel(size(opts.color)) == 2 && ...
        all(size(opts.color) == [2, 3]) && ...
       ~any(isnan(opts.color(:))) && ...
        all(opts.color(:) >= 0 & opts.color(:) <= 255)
        iCBG = opts.color(1, :);
        iCFG = opts.color(2, :);
        if any(opts.color(:) > 1)
            iCBG = iCBG ./ 255;
            iCFG = iCFG ./ 255;
        end
    else
        iCBG = [ 16,  32, 192] / 255;
        iCFG = [128, 192, 240] / 255;
    end

    % put together bar image
    imgpdata(numlines, 1, 3) = uint8(0);
    for colc = 1:numcols
        colb = (numcols - colc) / numcols;
        colf = colc / numcols;
        ridx = [colc (numlines + 1 - colc)];
        imgpdata(ridx, 1, 1) = ...
            uint8(fix(255 * (colb * iCBG(1) + colf * iCFG(1)) + 0.5));
        imgpdata(ridx, 1, 2) = ...
            uint8(fix(255 * (colb * iCBG(2) + colf * iCFG(2)) + 0.5));
        imgpdata(ridx, 1, 3) = ...
            uint8(fix(255 * (colb * iCBG(3) + colf * iCFG(3)) + 0.5));
    end

    % create first axes and add image
    AxProp = struct( ...
        'Parent',   hFig, ...
        'Color',    'none', ...
        'Units',    'normalized', ...
        'Position', [0.05, 0.15, 0.9, 0.7]);
    hBarAx = axes(AxProp);
    image(imgpdata, 'Parent', hBarAx);
    set(hBarAx, 'Visible', 'off');

    % create second axes for outline and caption
    hLabAx = axes(AxProp);
    if mlversion(1) < 8 || ...
        (mlversion(2) == 8 && ...
         mlversion(2) < 4)
        line(...
            [0 1 1 0 0], [0 0 1 1 0], ...
            'Parent',    hLabAx, ...
            'EraseMode', 'none', ...
            'Color',     iCBG * 0.75, ...
            'Visible',   'on');
    else
        line(...
            [0 1 1 0 0], [0 0 1 1 0], ...
            'Parent',    hLabAx, ...
            'Color',     iCBG * 0.75, ...
            'Visible',   'on');
    end

    CapProp = struct('Parent', hLabAx);
    CapProp.Units           = 'normalized';
    CapProp.Position        = [0.02 0.52];
    CapProp.HorizontalAlign = 'left';
    CapProp.VerticalAlign   = 'middle';
    CapProp.FontSize        = 11;
    CapProp.FontWeight      = 'demi';
    CapProp.FontUnits       = 'points';
    CapProp.Interpreter     = 'none';
    if prod(iCFG) < 0.125 && ...
        all(iCFG < 0.25)
        CapProp.Color       = [0.9325 0.9325 0.9325];
    else
        CapProp.Color       = iCBG * 0.125;
    end

    % generate text object
    hBarTxt = text(0.5, 0.5, iCap, CapProp);
    if numel(get(hBarTxt, 'Position')) == 3
        set(hBarTxt, 'Position', [0.02, 0.52, 0]);
    end

    % get size, set Units back to original units and init progress
    set(hLabAx, 'Units', 'normalized', 'Visible', 'off');

    % fill class object
    p = struct;
    p.hFig = hFig;
    p.hBar = hBarAx;
    p.hLab = hBarTxt;
    p.pBar = get(p.hBar, 'Position');
    p.sTitle = 'Progress...';
    p.vCur = 0;
    p.vDst = 1;
    p.vMax = 1;
    p.vMin = 0;
    p.vNow = 86400 * now;
    p.vPID = randn(1, 1);

    % create class object
    p = class(p, 'xprogress');

    % set bar width
    p = xprogress(p, p.vCur);

    % return
    return;

end

% empty call
if nargin == 1
    p = p.vCur;
    return;
end

% command call
if nargin > 2 && ...
    ischar(pval) && ...
    any(strcmpi(pval(:)', ...
        {'close', 'setcaption', ...
         'setposition', 'setsize', ...
         'settitle', 'setvisible'}))

    % what command ?
    switch (lower(pval(:)'))

        % close fig
        case {'close'}
            try
                xpconfig.lastpos = get(p.hFig, 'Position');
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
            end
            delete(p.hFig);

        % caption
        case {'setcaption'}
            try
                set(p.hLab, 'String', ptxt);
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
            end

        % position
        case {'setposition'}
            try
                if numel(ptxt) ~= 4
                    fpos = get(p.hFig, 'Position');
                    set(p.hFig, 'Position', [ptxt(1:2), fpos(3:4)]);
                else
                    set(p.hFig, 'Position', ptxt);
                end
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
            end

        % size
        case {'setsize'}
            try
                fpos = get(p.hFig, 'Position');
                if numel(ptxt) < 2
                    set(p.hFig, 'Position', [fpos(1:2), ptxt(1), fpos(4)]);
                else
                    set(p.hFig, 'Position', [fpos(1:2), ptxt(1:2)]);
                end
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
            end

        % title
        case {'settitle'}
            try
                p.sTitle = ptxt;
                set(p.hFig, 'Name', ptxt);
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
            end

        % visibility
        case {'setvisible'}
            try
                set(p.hFig, 'Visible', ptxt);
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
            end
    end

    % set back
    if ~isempty(inputname(1))
        try
            assignin('caller', inputname(1), p);
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
        end
    end

    % return
    drawnow;
    return;
end

% objects must still exist
if ~ishandle(p.hFig) || ~ishandle(p.hBar) || ~ishandle(p.hLab)
    warning('neuroelf:xprogress:badHandle', ...
        'Handle of figure, bar or text object disappeared.');
    return;
end

% genuine call with initial settings?
if nargin > 1 && isa(pval, 'double') && numel(pval) == 1 && ...
   ~isinf(pval) && ~isnan(pval)

    % make sure value is OK
    pval = min(p.vMax, max(p.vMin, pval));
    p.vCur = pval;

    % set size
    pval = max(p.vMin + eps, pval);
    pval = (pval - p.vMin) / p.vDst;
    etime = 86400 * now - p.vNow;
    ttime = ceil(etime / (pval + eps) + 0.1);
    etime = ceil(etime + 0.1);
    ttime = max(ttime, etime);
    tlstr = p.sTitle;
    if numel(tlstr) > 40
        tlstr = [tlstr(1:19) '...' tlstr(end-18:end)];
    end
    if pval > 0.001
        if ttime < 90
            ttstr = sprintf('%02d/%02ds', etime, ttime);
        elseif ttime < 3600
            ttstr = sprintf('%02d:%02d/%02d:%02d', ...
                floor(mod(etime, 3600) / 60), mod(etime, 60), ...
                floor(mod(ttime, 3600) / 60), mod(ttime, 60));
        else
            ttstr = sprintf('%02d:%02d:%02d/%02d:%02d:%02d', ...
                floor(etime / 3600), floor(mod(etime, 3600) / 60), mod(etime, 60), ...
                floor(ttime / 3600), floor(mod(ttime, 3600) / 60), mod(ttime, 60));
        end
    else
        ttstr = sprintf('%02d:%02d:%02d/??:??:??', ...
            floor(etime / 3600), floor(mod(etime, 3600) / 60), mod(etime, 60));
    end
    set(p.hBar, 'Position', [p.pBar(1:2), p.pBar(3) * pval, p.pBar(4)]);
    set(p.hFig, 'Name', sprintf('%.1f%% (%s) - %s', 100 * pval, ttstr, tlstr));
    % anything else
    if nargin > 2 && ischar(ptxt) && ~isempty(ptxt)
        set(p.hLab, 'String', ptxt(:)');
    end
    if nargin > 3 && ischar(pstate) && ...
        any(strcmpi(pstate(:)', {'close', 'hidden', 'visible'}))
        switch (lower(pstate(:)'))
            case 'close'
                xprogress(p, 'close', []);
                return;
            case 'hidden'
                set(p.hFig, 'Visible', 'off');
            case 'visible'
                set(p.hFig, 'Visible', 'on');
        end
    end
    if nargin > 5 && ...
        isa(pmin, 'double') && ...
        numel(pmin) == 1 && ...
        ~isinf(pmin) && ...
        ~isnan(pmin) && ...
        isa(pmax, 'double') && ...
        numel(pmax) == 1 && ...
        ~isinf(pmax) && ...
        ~isnan(pmax) && ...
        pmax > pmin
        p.vMin = pmin;
        p.vMax = pmax;
        p.vDst = pmax - pmin;
        p.vCur = min(pmax, max(pmin, p.vCur));
        xprogress(p, p.vCur);
    end

    % update screen
    drawnow;

    % set back
    if ~isempty(inputname(1))
        try
            assignin('caller', inputname(1), p);
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
        end
    end
end
