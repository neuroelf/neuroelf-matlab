function o = progress(xo, iStr, varargin)
%XFIGURE::PROGRESS  Update and return progress bar control.
%   PROGRESS(PBAR, PROGVALUE) updates the progress bar PBAR with value
%   PROGVALUE.
%
%   P = PROGRESS(PBAR) returns the current progress value.

% Version:  v1.1
% Build:    16041117
% Date:     Apr-11 2016, 5:03 PM EST
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

% single xprogress argument
if numel(xo) ~= 1 || xo.T ~= 2 || ~strcmpi(xo.X.loadprops.xtype, 'xprogress')
    error('neuroelf:xfigure:invalidObjectType', ...
        'Progress only valid for UIControl/XProgress objects.');
end

% empty progress
if nargin < 2 || ~isnumeric(iStr) || isempty(iStr)
    refresh(get(xo.H, 'Parent'));
    o = xo.X.uicprops.progress;
    return;
end

% set width?
if nargin > 2 && ischar(varargin{1}) && any(strcmpi(varargin{1}(:)', {'newwidth', 'newheight'}))
    % po = xfigure(get(xo.H, 'Parent'));
    if strcmpi(varargin{1}(:)', 'newwidth')
        xo.X.loadprops.Position(3) = iStr;
        pbpos = get(xo.H, 'Position');
        set(xo.H, 'Position', [pbpos(1:2), iStr, pbpos(4)]);
    elseif strcmpi(varargin{1}(:)', 'newheight')
        xo.X.loadprops.Position(4) = iStr;
        pbpos = get(xo.H, 'Position');
        set(xo.H, 'Position', [pbpos(1:3), iStr]);
    end
    if nargin < 4
        o = progress(xo, xo.X.uicprops.progress);
    else
        o = progress(xo, xo.X.uicprops.progress, [], varargin{2:end});
    end
    return;
end

% get handle
hPatch = xo.X.uicprops.xchildren(2);

% calculate progress position
if isinf(iStr(1)) || isnan(iStr(1))
    PrgPos = xo.X.uicprops.progress;
else
    PrgPos = max(0, min(1, (iStr(1) - xo.X.loadprops.MinMaxTop(1)) / ...
             (xo.X.loadprops.MinMaxTop(2) - xo.X.loadprops.MinMaxTop(1))));
    xo.X.uicprops.progress = PrgPos;
end

% return value
o = xo.X.uicprops.progress;

% update?
if xo.X.uicprops.nextupdate > now
    return;
end
xo.X.uicprops.nextupdate = now + xo.X.uicprops.updrate;

% round progress bars:
% we must get the position from the main axes object and then
% use this to calculate the image axes' position
switch (lower(xo.X.loadprops.ProgType))
    case 'round'
        try
            hBarAx = get(xo.X.uicprops.xchildren(1), 'Parent');
            set(hBarAx, 'Units', 'pixels');
            oPos = xo.X.loadprops.Position;
            if strcmp(xo.X.loadprops.ProgDir, 'y')
                oPos(4) = oPos(4) * PrgPos + eps;
            else
                oPos(3) = oPos(3) * PrgPos + eps;
            end
            set(hBarAx, 'Position', oPos);
        catch xfigerror
            warning(xfigerror.message);
        end

    % flat progress bars:
    case 'flat'
        try
            % horiz or vert progress bar
            if strcmp(xo.X.loadprops.ProgDir, 'y')
                set(hPatch, 'Position', [0, 0, 1, min(1, (PrgPos + eps))]);
            else
                set(hPatch, 'Position', [0, 0, min(1, (PrgPos + eps)), 1]);
            end
        catch xfigerror
            warning(xfigerror.message);
        end

    % otherwise, bail out
    otherwise
        error('neuroelf:xfigure:internalError', ...
            'Invalid ProgressBar type: %s.', xo.X.loadprops.ProgType);
end

% labeling
newlabel = {};
if nargin > 2 && ischar(varargin{1}) && numel(xo.X.uicprops.xchildren) > 2
    newlabel = {varargin{1}(:)'};
    try
        set(xo.X.uicprops.xchildren(end), 'String', newlabel{1});
    catch xfigerror
        warning(xfigerror.message);
    end
end
if nargin < 4 || ~islogical(varargin{2}) || numel(varargin{2}) ~= 1 || varargin{2}
    drawnow;
end

% additional updating?
upbar = get(xo.H, 'UserData');
if iscell(upbar) && ~isempty(upbar) && numel(upbar{1}) == 1 && isa(upbar{1}, 'function_handle')
    feval(upbar{:}, PrgPos, newlabel{:});
end
