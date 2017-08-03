function pgnum = showpage(xo, iStr, varargin)
%XFIGURE::SHOWPAGE  Show page of a figure (internal set of VGroup elements).
%   SHOWPAGE(FIG, P) shows page P of figure FIG.

% Version:  v1.1
% Build:    16041116
% Date:     Apr-11 2016, 4:40 PM EST
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

% for figures
if numel(xo) == 1 && xo.T == 1

    % character page spec
    if ischar(iStr) && ~isempty(iStr) && ~isempty(xo.X.figprops.pages)

        % for specific tokens
        switch lower(iStr(:)')
            case 'cur'
                pgnum = xo.X.figprops.cpage;
                return;
            case 'max'
                pgnum = max(xo.X.figprops.pages);
                return;
            case 'min'
                pgnum = min(xo.X.figprops.pages);
                return;
        end

        % parse page
        try
            tpage = str2double(iStr);
            if isempty(tpage) || isnan(tpage(1)) || isinf(tpage(1))
                error('INVALID_PAGE_NUMBER_STRING');
            end
        catch xfigerror
            neuroelf_lasterr(xfigerror);
            error('neuroelf:xfigure:badArgument', 'Invalid page (string) requested.');
        end
        if ~any(iStr(1) == '0123456789')
            tpage = max(xo.X.figprops.cpage, 1) + tpage(1);
        end
        iStr = tpage(1);
    elseif ~isnumeric(iStr) || isempty(iStr) || isnan(iStr(1)) || isinf(iStr(1))
        error('neuroelf:xfigure:badArgument', ...
            'Invalid page (datatype/content) requested.');
    end
    pgnum  = fix(max(1, min(max(xo.X.figprops.pages), iStr(1))));
    xo.X.figprops.cpage = pgnum;
    callhdl = xo;

% for uicontrols
elseif numel(xo) == 1 && xo.T == 2
    pgnum = get(xo.H, 'Value');
    iFObj = xfigure(get(xo.H, 'Parent'));
    pgspec = iFObj.X.figprops.pages;
    if pgnum > length(pgspec)
        warning('neuroelf:xfigure:badShowPageControl', ...
            'Invalid page selected -> wrong dropdown?');
        return;
    end
    if pgnum > 0
        pgnum = iFObj.X.figprops.pages(pgnum);
    end
    iFObj.X.figprops.cpage = pgnum;
    callhdl = iFObj;

% otherwise
else
    error('neuroelf:xfigure:invalidObjectType', ...
        'ShowPage is only valid for figures and uicontrols.');
end

% make calls to SetGroupVisible
if nargin > 2 && ischar(varargin{1}) && ~strcmpi(varargin{1}, 'refresh')
    refreshwin = {};
else
    refreshwin = {'refresh'};
end
setgroupvisible(callhdl, 'UICPage_any', 'off');
setgroupvisible(callhdl, 'UICPage_all');
setgroupvisible(callhdl, ['UICPage' num2str(pgnum)], 'on', refreshwin{:});
