function [tdreply, st] = tdclient(varargin)
% tdclient  - send request to internet TDdatabase / local TD
%
%       FORMAT:     tdreply = tdclient( x,y,z  [, kind [, cubesize]])
%       alt.        tdreply = tdclient([x,y,z] [, kind [, cubesize]])
%
%       FORMAT:     tdclient('format', formstring, formsep [, extraline])
%
% Input fields:
%       x           x coordinate of request
%       y           y coordinate of request
%       z           z coordinate of request
%       kind        string, {'TalLabel'}, 'Cube', 'NGM', 'RNGM'
%       cubesize    size of voxelcube to search, {7} (from 3 through 11)
%
%       formstring  specifies the output format for RNGM call
%       formsep     sets the field seperator
%       extraline   enable/disable extra line feed in output
%
%   NOTE:           in the alternative calling form, a ?x3 matrix can be
%                   processed, then a cell array will be returned!
%
% when calling with kind='Cube' or kind='NGM' a cubesize of 7mm is
% assumed if not specified in the argument list; in NGM mode, only those
% points will be returned that lie within a radius of cubesize !
%
% as a special feature, NGM offers a 1x2 cubesize argument, where the
% second value must be greater than 1 to show the first n occurences of
% the search, in command mode, seperate with a '#' (hash)
%
% tdclient can be called as a MATLAB command with strings as input!
% here are some examples for this kind of request:
%
% tdclient -21 -24  12          returns what lies at -21,-24,12
% tdclient   5  17  -5 Cube 9   returns the number of voxels for each label
%                                 found in a 9x9x9 cube around 5,17,-5
% tdclient -22 -21  22 NGM 5    finds the nearest Gray Matter from
%                                 -22,-21,22; with a radius of 5mm
% tdclient -22 -21  22 NGM 7#5  find the 5 nearest distinct Gray Matter
%                                 occurences from the same point, radius 7mm
%
% when run in command mode, the return value will always be formatted
%
% due to the fact that extraline is a numeric value, the formatting call
% must be made as a function:
%
% tdclient('format','$xc;$yc;$zc;$diff;$xd;$yd;$zd;$t1;$t3;$t5;',';',0);
%
% where $xc, $yc, $zc are replaced with the actual coordinates of the
% gray matter found, $diff is the calculated distance between the
% requested and the actual coordinates, and $xd, $yd, $zd are the
% differences for each coordinate. $t1 through $t5 are the text patterns
% for each group of location descriptors.
%
% See also tdlocal2.

% Version:  v0.9d
% Build:    14060710
% Date:     Jun-07 2014, 10:25 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2014, Jochen Weber
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

% persistent variable
persistent td_iset;
if isempty(td_iset)

    % test required local version of TD database
    try
        tdlocal2(2, 0, 0, 0);
    catch ne_eo;
        error( ...
            'neuroelf:TDError', ...
            'Error using tdlocal2: ''%s''', ...
            ne_eo.message ...
        );
    end
    td_iset.gui = [];
end

% enough arguments ?
if nargin < 1
    error( ...
        'neuroelf:TooFewArguments',...
        'Too few arguments.' ...
    );
end

% preset output
if nargout > 0
    tdreply = '';
end

% settings
tddecimal = '.';
tdformat  = '$xc;$yc;$zc;$diff;$xd;$yd;$zd;$t1;$t2;$t3;$t4;$t5;';
tdforsep  = ';';
tdforxln  = false;
indent    = '                 ';
lfc       = char(10);

% - set initial values, settings...
kinds  = struct(   ...
    'tallabel', 2, ...
    'cube',     3, ...
    'ngm',      5, ...
    'rngm',    -5  ...
);

if ischar(varargin{1}) && ...
    strcmpi('visual', varargin{1})
    if nargin == 1 && ...
        isempty(td_iset.gui)
        try
            hWnd = xfigure([neuroelf_path('tfg') '/tdclient.tfg']);
            hWnd.HandleVisibility = 'callback';
            td_iset.gui.hWnd = hWnd;
            td_iset.gui.hnds = hWnd.TagStruct;
            tdclient('visual', 'clear');
        catch ne_eo;
            warning( ...
                'neuroelf:xfigureError', ...
                'Couldn''t open GUI: %s', ...
                ne_eo.message ...
            );
            return;
        end
    elseif nargin == 1
        return;
    elseif ischar(varargin{2})
        switch (lower(varargin{2}))
            case {'clear'}
                try
                    td_iset.gui.hnds.ED_tdclient_input.String = '';
                    td_iset.gui.hnds.ED_tdclient_result.String = '';
                    td_iset.gui.hnds.RB_tdclient_TALLabel.DoCallback;
                    td_iset.gui.hnds.ED_tdclient_SrcRange.String = '7';
                    td_iset.gui.hnds.ED_tdclient_NGMoccur.String = '5';
                catch ne_eo;
                    neuroelf_lasterr(ne_eo);
                end
            case {'close'}
                try
                    td_iset.gui = [];
                    hWnd = xfigure('Wnd_tdclient');
                    hWnd.Delete;
                catch ne_eo;
                    neuroelf_lasterr(ne_eo);
                end
            case {'send'}
                try
                    try
                        input = evalin('base', ...
                            ['[' td_iset.gui.hnds.ED_tdclient_input.String ']']);
                    catch ne_eo;
                        uiwait(warndlg(sprintf('Error parsing input: %s', ne_eo.message), ...
                            'tdclient - error', 'modal'));
                        return;
                    end
                    if ~any(size(input) == 3)
                        uiwait(warndlg('Input must be Cx3 or 3xC.', ...
                            'tdclient - error', 'modal'));
                        return;
                    end
                    if td_iset.gui.hnds.RB_tdclient_TALLabel.IsActive
                        retval = tdclient(input);
                    elseif td_iset.gui.hnds.RB_tdclient_CubeSearch.IsActive
                        try
                            srcrange = round(str2double( ...
                                td_iset.gui.hnds.ED_tdclient_SrcRange.String));
                        catch ne_eo;
                            neuroelf_lasterr(ne_eo);
                            srcrange = 7;
                        end
                        retval = tdclient(input, 'Cube', srcrange);
                    else
                        try
                            srcrange = round(str2double( ...
                                td_iset.gui.hnds.ED_tdclient_SrcRange.String));
                            numhits = round(str2double( ...
                                td_iset.gui.hnds.ED_tdclient_NGMoccur.String));
                        catch ne_eo;
                            neuroelf_lasterr(ne_eo);
                            srcrange = 7;
                            numhits = 5;
                        end
                        retval = tdclient(input, 'NGM', [srcrange, numhits]);
                    end
                    if iscell(retval)
                        retval = gluetostring(retval);
                    end
                    td_iset.gui.hnds.ED_tdclient_result.String = retval;
                catch ne_eo;
                    neuroelf_lasterr(ne_eo);
                end
        end
    end
	return;
end

% argument check
if ischar(varargin{1}) && ...
    any(varargin{1}(:) == ',')
    st = 0;
    tdargs = splittocell(varargin{1}, ',', 1);
    if any(tdargs{1} == ':')
        tdargsx = splittocell(tdargs{1}, ':', 1);
        tdargs{1} = tdargsx{1};
        tdargs{end+1} = tdargsx{2};
    end
    tdreply = tdlocal2(tdargs{1:end});
    tdreply = [varargin{1} ': ' tdreply];

else
    if ischar(varargin{1})
        if nargin < 3
            error( ...
                'neuroelf:TooFewArguments',...
                'Too few arguments. Try ''help %s''.',...
                mfilename ...
            );
        end
        varargin{1} = str2double(varargin{1});
        if isempty(varargin{1})
            error( ...
                'neuroelf:BadArgument',...
                'Bad argument.' ...
            );
        end
        if ischar(varargin{2}), varargin{2} = str2double(varargin{2}); end
        if ischar(varargin{3}), varargin{3} = str2double(varargin{3}); end
        if nargin > 4 && ...
            ischar(varargin{5})
            eval(['varargin{5} = [' strrep(varargin{5},'#',',') '];']);
        end
        charout = 1;
    else
        charout = 0;
    end
    fromspm = 0;
    if isstruct(varargin{1}) && ...
        isfield(varargin{1}, 'dat')
        fromspm = 1;
        tdata   = varargin{1}.dat;
        if ~iscell(tdata)
            error( ...
                'neuroelf:BadArgument',...
                'Bad argument.' ...
            );
        end
        tdat  = tdata(:,size(tdata,2));
        tdats = size(tdat,1);
        tci   = zeros(tdats,3);
        tke   = cell(1,tdats);
        tze   = cell(1,tdats);
        for tc = 1:tdats
            [tci(tc,1:3)] = round([tdat{tc}(1), tdat{tc}(2), tdat{tc}(3)]);
            tke{tc} = sprintf('%.0f', tdata{tc,4});
            tze{tc} = strrep(sprintf('%.4f', tdata{tc,9}), '.', tddecimal);
        end
        varargin{1} = tci;
    end
    if ~isnumeric(varargin{1})
        error( ...
            'neuroelf:BadArgument',...
            'Bad argument.' ...
        );
    end
    if length(varargin{1}) < 3
        coords = [varargin{1:3}];
        if nargin > 3
            ikind = varargin{4};
        else
            ikind = 2;
        end
        if nargin > 4
            isize = varargin{5};
        else
            isize = 7;
        end
    else
        coords = varargin{1};
        if nargin > 1
            ikind = varargin{2};
        else
            ikind = 2;
        end
        if nargin > 2
            isize = varargin{3};
        else
            isize = 7;
        end
    end
    if ischar(ikind)
        if isfield(kinds,lower(ikind))
            ikind = kinds.(lower(ikind));
        else
            ikind = 2;
        end
    end
    if ischar(isize)
        try
            isize = eval(['[' strrep(isize, '#', ',') ']']);
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            isize = [7, 5];
        end
    end
    if (ikind < 2 || ikind > 5) && ...
        ikind ~= -5
        ikind = 2;
    end
    if ~isnumeric(isize)
        isize = 7;
    end
    if mod(isize(1),2) ~= 1 && ...
        ikind ~= 5 && ...
        ikind ~= -5
        isize(1) = isize(1) + 1;
    end
    if isize < 3
        isize = 3;
    end
    if isize > 11
        isize = 11;
    end

    numcoords = size(coords, 1);
    tdreplya = cell(1, numcoords);
    for rc = 1:numcoords
        if ikind > 0 && ikind < 3
            command = sprintf('%d,%d,%d,%d',     ikind(1), coords(rc,1:3));
        elseif ikind == 3 || ((ikind == 5 || ikind == -1) && isize(1) == fix(isize(1)))
            command = sprintf('3:%d,%d,%d,%d',   isize(1), coords(rc,1:3));
        elseif (ikind == 5 || ikind == -5)
            command = sprintf('5:%.2f,%d,%d,%d', isize(1), coords(rc,1:3));
        end
        if ikind == 2
            st = 0;
            tdreply = tdlocal2(ikind(1), coords(rc, 1), coords(rc, 2), coords(rc, 3));
        else
            st = 0;
            tdreply = tdlocal2(ikind(1), coords(rc, 1), coords(rc, 2), coords(rc, 3), isize);
        end
        if ikind == 3 || ...
            ikind == 5
            tdreplyc = splittocell(tdreply,':',0);
            tdreply  = '';
            for cc = 2:2:length(tdreplyc)
                tdreplys = [tdreplyc{cc-1} ': ' tdreplyc{cc}];
                myindent = indent;
                if cc == 2
                    myindent(1:length(command)) = command(1:end);
                end
                tdreply = [tdreply myindent tdreplys lfc];
            end
            tdreplya{rc} = tdreply;
        elseif ikind == -5
            if ~isempty(tdforsep)
                tdforsep = tdforsep(1);
            else
                tdforsep = ';';
            end
            lformat(1:length(find(tdformat == tdforsep))) = tdforsep;
            if fromspm
                sformat(1:7) = tdforsep;
                tdreplyr = sprintf('%d%s%d%s%d%s%s%s%s%s%s%s%s%s%s%s%s', ...
                    coords(rc,1), tdforsep, ...
                    coords(rc,2), tdforsep, ...
                    coords(rc,3), tdforsep, ...
                    'rngm', tdforsep, ...
                    sprintf('%d:', isize), tdforsep, ...
                    tke{rc}, tdforsep, ...
                    tze{rc}, tdforsep, ...
                    lformat, lfc);
            else
                sformat(1:5) = tdforsep;
                tdreplyr = sprintf('%d%s%d%s%d%s%s%s%s%s%s%s%s', ...
                    coords(rc,1), tdforsep, ...
                    coords(rc,2), tdforsep, ...
                    coords(rc,3), tdforsep, ...
                    'rngm',       tdforsep, ...
                    sprintf('%d:', isize), tdforsep, ...
                    lformat, lfc);
            end
            if iscell(tdreply)
                for cc = 1:length(tdreply)
                    myformat = tdformat;
                    myformat = strrep(myformat, '$xc', ...
                                      sprintf('%d',tdreply{cc}{1}(1)));
                    myformat = strrep(myformat, '$yc', ...
                                      sprintf('%d',tdreply{cc}{1}(2)));
                    myformat = strrep(myformat, '$zc', ...
                                      sprintf('%d',tdreply{cc}{1}(3)));
                    myformat = strrep(myformat, '$xd', ...
                                      sprintf('%d',tdreply{cc}{2}(1)));
                    myformat = strrep(myformat, '$yd', ...
                                      sprintf('%d',tdreply{cc}{2}(2)));
                    myformat = strrep(myformat, '$zd', ...
                                      sprintf('%d',tdreply{cc}{2}(3)));
                    myformat = strrep(myformat, '$diff', ...
                                      strrep(sprintf('%d', ...
                                          tdreply{cc}{3}(1)),'.',tddecimal));
                    myformat = strrep(myformat, '$t1', ...
                                      tdreply{cc}{4}{1});
                    myformat = strrep(myformat, '$t2', ...
                                      tdreply{cc}{4}{2});
                    myformat = strrep(myformat, '$t3', ...
                                      tdreply{cc}{4}{3});
                    myformat = strrep(myformat, '$t4', ...
                                      tdreply{cc}{4}{4});
                    myformat = strrep(myformat, '$t5', ...
                                      tdreply{cc}{4}{5});
                    tdreplyr = [tdreplyr sformat myformat lfc];
                end
            end
            if all(tdforxln == 1)
                tdreplya{rc} = [tdreplyr sformat lformat lfc];
            else
                tdreplya{rc} = tdreplyr(1:(end-length(lfc)));
            end
        else
            command  = [command ':'];
            myindent = indent;
            myindent(1:length(command)) = command;
            tdreplya{rc} = [myindent tdreply];
        end
    end

    if rc == 1 && fromspm == 0
        tdreply = tdreplya{1};
    else
        tdreply = tdreplya;
    end
    if charout, tdreply = char(tdreply); end
end
