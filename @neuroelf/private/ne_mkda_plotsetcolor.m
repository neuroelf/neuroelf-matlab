% PUBLIC FUNCTION ne_mkda_plotsetcolor: set color of points
function varargout = ne_mkda_plotsetcolor(varargin)

% Version:  v0.9c
% Build:    13012611
% Date:     Mar-07 2012, 12:01 PM EST
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
ccfg = ne_gcfg.fcfg;
ch = ne_gcfg.h.MKDA.h;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% get content of PLP dropdown
plps = ch.PLPs;
plpud = plps.UserData;
plpid = plps.Value;
try
    plp = plpud{plpid, 3};
    if numel(plp) ~= 1 || ...
       ~isxff(plp, 'plp')
        error( ...
            'neuroelf:GUI:BadPLP', ...
            'Bad PLP object.' ...
        );
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    return;
end

% no Color column
cn = plp.ColumnNames;
cn = cn(:);
if ~any(strcmpi(cn, 'color'))
    plp.ColumnNames{end+1} = 'Color';
    plp.NrOfColumns = numel(plp.ColumnNames);
    plp.Points(:, plp.NrOfColumns) = 1;
    plp.RunTimeVars.ColumnIsText.Color = false;
    cn = plp.ColumnNames;
    cn = cn(:);
end
if plp.NrOfColors < 1
    plp.NrOfColors = 1;
    plp.Colors = [255, 0, 0];
end
ccol = findfirst(strcmpi(cn, 'color'));

% which column
if nargin < 3 || ...
   ~ischar(varargin{3}) || ...
    isempty(varargin{3})
    columnid = [];
else
    columnid = findfirst(strcmpi(cn, varargin{3}(:)'));
end
if isempty(columnid)
    columnid = listdlg( ...
        'ListString',    [cn; {'SingleColor'}], ...
        'SelectionMode', 'single', ...
        'ListSize',      [max(320, min(640, 10 * size(char(cn), 2))), 320], ...
        'InitialValue',  1, ...
        'Name',          'NeuroElf - user input', ...
        'PromptString',  'Please select a column to set a color code...');
    if numel(columnid) ~= 1 || ...
        columnid < 1 || ...
        columnid > (numel(cn) + 1)
        return;
    end
end

% single color
if columnid > numel(cn)
    try
        plp.Colors(1, :) = colorpicker(plp.Colors(:, 1), {'PLP points color'});
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
    return;
end

% get RunTimeVars
rtv = plp.RunTimeVars;

% and unique values
cdat = plp.Points(:, columnid);
uval = unique(cdat);

% if the column is a text column
cit = rtv.ColumnIsText.(cn{columnid});
if cit

    % we need the labels
    labels = plp.Labels(:);

    % get unique labels
    ulab = labels(uval);

    % get colors associated with labels
    if size(plp.Colors, 1) < uval(end)
        plp.Colors(end+1:uval(end), :) = ...
            floor(256 .* rand(uval(end) - size(plp.Colors, 1), 3));
        plp.NrOfColors = size(plp.Colors, 1);
    end
    ucol = plp.Colors(uval, :);

% column isn't text
else

    % create labels
    ulab = cell(numel(uval), 1);
    for lc = 1:numel(ulab)
        ulab{lc} = sprintf('%g', uval(lc));
    end

    % and generate colors
    ucol = floor(256 .* rand(numel(ulab), 3));
end

% request colors
ncol = colorpicker(ucol, ulab);

% nothing changed
if all(ncol(:) == ucol(:)) && ...
   ~cit

    % don't do anything
    return;
end

% column is text?
if cit

    % set into correct places
    plp.Colors(uval, :) = ncol;

    % and then set values
    plp.Points(:, ccol) = cdat;

% no text
else

    % overwrite the first N values
    plp.Colors(1:size(ncol, 1), :) = ncol;
    plp.NrOfColors = size(plp.Colors, 1);

    % then set values 1 through N to color column
    for lc = 1:numel(uval)
        plp.Points(cdat == uval(lc), ccol) = lc;
    end
end

% update saved flag
plp.RunTimeVars.Saved = false;
