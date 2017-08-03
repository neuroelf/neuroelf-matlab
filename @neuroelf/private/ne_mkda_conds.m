% PUBLIC FUNCTION ne_mkda_conds: deal with conditional statements
function varargout = ne_mkda_conds(varargin)

% Version:  v0.9c
% Build:    11120120
% Date:     Nov-08 2011, 1:08 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, Jochen Weber
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
ch = ne_gcfg.h.MKDA.h;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% make sure that 3rd argument is valid action
if nargin < 3 || ...
   ~ischar(varargin{3}) || ...
    isempty(varargin{3}) || ...
   ~any(strcmp(varargin{3}(:)', {'add', 'addpar', 'del', 'delpar', 'sel', 'set'}))
    return;
end
action = varargin{3}(:)';

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

% get current condition list and index
condparts = ch.CndParts.String;
if ~iscell(condparts)
    condparts = cellstr(condparts);
end
condidx = ch.CndParts.Value;

% depending on action
switch (action)

    % adding a new particle (or re-setting current one)
    case {'add', 'set'}

        % for setting, we need single selection
        if strcmp(action, 'set') && ...
            numel(condidx) ~= 1
            return;
        end

        % get linkage, column namem operator and operand
        condand = (ch.CndAnd.Value > 0);
        condcol = ch.CndColumn.String{ch.CndColumn.Value};
        condoper = ch.CndOperator.Value;
        condoprnd = ch.CndFlexiOp.String;

        % textcolumn
        if ch.CndColumn.UserData

            % equality
            if mod(condoper, 2) == 1
                condop = '==';
            else
                condop = '~=';
            end

            % wildcarding
            switch condoper

                % equality
                case {1, 2}
                    condoprnd = [  '''' condoprnd ''''  ];

                % contains
                case {3, 4}
                    condoprnd = ['''.*' condoprnd '.*'''];

                % begins
                case {5, 6}
                    condoprnd = [  '''' condoprnd '.*'''];

                % ends
                case {7, 8}
                    condoprnd = ['.*''' condoprnd ''''  ];
            end

        else

            % operator
            switch condoper

                % equals
                case {1}
                    condop = '==';

                % doesn't equal
                case {2}
                    condop = '~=';

                % less than
                case {3}
                    condop = '<';

                % less than or equal
                case {4}
                    condop = '<=';

                % greater than
                case {5}
                    condop = '>';

                % greater than or equal
                case {6}
                    condop = '>=';
            end
        end

        % entire condition
        if isempty(condparts) || ...
           (strcmp(action, 'set') && ...
            any(condidx == 1))
            cond = sprintf('   $%s %s %s ', condcol, condop, condoprnd);
        elseif condand
            cond = sprintf(' & $%s %s %s ', condcol, condop, condoprnd);
        else
            cond = sprintf(' | $%s %s %s ', condcol, condop, condoprnd);
        end

        % add to list
        if strcmp(action, 'add')
            condparts{end+1} = cond;

        % set particle
        else

            % where in particle
            ocond = deblank(condparts{condidx});
            cfidx = findfirst(ocond == '$');
            clidx = numel(ocond);
            while ocond(clidx) == ')'
                clidx = clidx - 1;
            end
            condparts{condidx} = [cond(1:3) ocond(4:(cfidx-1)) ...
                cond(4:end-1) ocond((clidx+1):numel(ocond))];
        end
        ch.CndParts.String = condparts;

    % remove particle(s)
    case {'del'}

        % only valid if index not empty
        if isempty(condidx)
            return;
        end

        % for a pre-test, get conditions that are to be combined
        cond = gluetostringc(condparts(condidx), ' ');

        % only valid if number of closing parentheses remains smaller!
        pdiff = cumsum(double(cond(:) == '(')) - cumsum(double(cond(:) == ')'));
        if any(pdiff < 0) || ...
            pdiff(end) ~= 0
            uiwait(warndlg('Cannot remove these particles due to parentheses.', ...
                'NeuroElf - error', 'modal'));
            return;
        end

        % remove particles
        condparts(condidx) = [];

        % if first particle was removed
        if any(condidx == 1) && ...
           ~isempty(condparts)

            % remove operator in front of now first particle
            condparts{1}(1:3) = '   ';
        end

        % then update
        ch.CndParts.Value = [];
        ch.CndParts.String = condparts;

    % add parantheses
    case {'addpar'}

        % only valid with index
        if numel(condidx) < 2
            return;
        end

        % for a pre-test, get conditions that are to be combined
        cond = gluetostringc(condparts(condidx), ' ');

        % only valid if number of closing parentheses remains smaller!
        pdiff = cumsum(double(cond(:) == '(')) - cumsum(double(cond(:) == ')'));
        if any(pdiff < 0)
            uiwait(warndlg('Cannot add parenthesis around these particles.', ...
                'NeuroElf - error', 'modal'));
            return;
        end

        % add parentheses
        condparts{condidx(1)} = ...
            [condparts{condidx(1)}(1:3) '(' condparts{condidx(1)}(4:end)];
        for cc = 2:numel(condidx)
            condparts{condidx(cc)} = ...
                [condparts{condidx(cc)}(1:3) ' ' condparts{condidx(cc)}(4:end)];
        end
        condparts{condidx(end)} = [deblank(condparts{condidx(end)}) ')'];

        % update
        ch.CndParts.String = condparts;

    % remove parantheses
    case {'delpar'}

        % only valid if index not empty
        if isempty(condidx)
            return;
        end

        % for a pre-test, get conditions that are to be combined
        cond = gluetostringc(condparts(condidx), ' ');

        % remove operator
        cond(1:3) = [];

        % find first space
        fspace = findfirst(cond == ' ');

        % only valid if number of closing parentheses remains smaller!
        pdiff = cumsum(double(cond(:) == '(')) - cumsum(double(cond(:) == ')'));
        if any(pdiff < 0) || ...
            pdiff(1) ~= 1 || ...
            pdiff(end-1) <= pdiff(end) || ...
            pdiff(end) > pdiff(fspace)
            uiwait(warndlg('Cannot remove any parentheses in this case.', ...
                'NeuroElf - error', 'modal'));
            return;
        end

        % remove parenthesis
        condparts{condidx(1)}(4) = [];
        for cc = 2:numel(condidx)
            if condparts{condidx(cc)}(4) == ' '
                condparts{condidx(cc)}(4) = [];
            end
        end
        condparts{condidx(end)}(end) = [];

        % update
        ch.CndParts.String = condparts;

    % select
    case {'sel'}

        % only valid for single selection
        if numel(condidx) ~= 1
            return;
        end

        % get selected particle
        cond = deblank(condparts{condidx});

        % remove unwanted parts
        cond(1:(findfirst(cond == '$')-1)) = [];
        while ~isempty(cond) && ...
            cond(end) == ')'
            cond(end) = [];
        end

        % nothing to do (error!)
        if isempty(cond)
            return;
        end

        % get column name
        [colname, cond] = strtok(cond(2:end), ' ');
        cond = ddeblank(cond);

        % try to select column
        colnames = ch.CndColumn.String;
        if ~iscell(colnames)
            colnames = cellstr(colnames);
        end
        colmatch = find(strcmpi(colname, colnames));
        if numel(colmatch) ~= 1
            return;
        end
        ch.CndColumn.Value = colmatch;
        ne_mkda_setcondcol;

        % update operator value as well
        [condop, condval] = strtok(cond, ' ');
        condval = ddeblank(condval);

        % for text columns
        if ch.CndColumn.UserData

            % operator
            if condop(1) == '='
                condoper = 1;
            else
                condoper = 2;
            end
            if strcmp(condval(2:3), '.*') && ...
                numel(condval) >= 4
                if strcmp(condval(end-2:end-1), '.*') && ...
                    numel(condval) >= 6
                    condval = condval(4:end-3);
                    condoper = condoper + 2;
                else
                    condval = condval(4:end-1);
                    condoper = condoper + 6;
                end
            else
                if strcmp(condval(end-2:end-1), '.*') && ...
                    numel(condval) >= 4
                    condval = condval(2:end-3);
                    condoper = condoper + 4;
                else
                    condval = condval(2:end-1);
                end
            end

        % for numeric columns
        else

            % operator
            switch (condop)
                case {'=='}
                    condoper = 1;
                case {'~='}
                    condoper = 2;
                case {'<'}
                    condoper = 3;
                case {'<='}
                    condoper = 4;
                case {'>'}
                    condoper = 5;
                case {'>='}
                    condoper = 6;
            end

            % sanitize condval
            if any(condval == ' ')
                condval = strtok(condval, ' ');
            end
        end
        ch.CndOperator.Value = condoper;
        ch.CndFlexiOp.String = condval;

        % try to match condval to listed values
        condvals = ch.CndStaticOp.String;
        if ~iscell(condvals)
            condvals = cellstr(condvals);
        end
        condvm = find(strcmpi(condval, condvals));
        if numel(condvm) ~= 1
            return;
        end
        ch.CndStaticOp.Value = condvm;

        % return early (don't update list!)
        return;

    % set (single particle)
    case {'set'}

        % only valid for single selection
        if numel(condidx) ~= 1
            return;
        end

        % get condition particle
        cond = condparts{condidx};

end

% update list
ne_mkda_listpoints(0, 0, 'upana');
