%ANY2ASCII  Pack simple as well as complex variables into ASCII string.
%   ASCII = ANY2ASCII(ANYVAR) converts the contents of the variable ANYVAR
%   into an ASCII representation, which can be used to either store the
%   content into a text file or for later use with EVAL.
%
%   The output will always be contained in a pair of [ and ] brackets,
%   whether or not ANYVAR is an array of a different size than 1x1.
%
%   ASCII = ANY2ASCII(ANYVAR, PRECISION) sets the precision for SPRINTF
%   for regular numeric (double) values. The default PRECISION is 8.
%
%   If precision is set to 'exact', ANY2ASCII uses the functions HXDOUBLE
%   and HXSINGLE to represent double-precision and single-precision float
%   values as the hexadecimal byte representation.
%
%   ANY2ASCII handles special values, such as (+/-)Inf, NaN, [], as well
%   as complex datatypes, including CELL and STRUCT arrays. In addition,
%   it also creates correctly sized empty arrays (e.g. 0x3x2 double), and
%   if the input is sparse, this will be reflected in the output as well.
%
%   For complex datatypes, the function is called recursively.
%
%   Importantly, if you wish a user-defined class to be supported, make
%   sure to overload the method (with supporting up to three arguments,
%   the third only being used for internal purposes).
%
%   See also DISP, EVAL, HXDOUBLE, HXSINGLE.

% Version:  v1.1
% Build:    16031110
% Date:     Mar-11 2016, 10:06 AM EST
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

% function definition
function asciirep = any2ascii(anyvar, varargin)

% persistent config
persistent a2acfg;
if isempty(a2acfg)
    a2acfg.classes = struct( ...
        'cell',    'l', ...
        'char',    'c', ...
        'double',  'd', ...
        'single',  'i', ...
        'int8',    '7', ...
        'int16',   '5', ...
        'int32',   '1', ...
        'int64',   '3', ...
        'logical', '9', ...
        'struct',  's', ...
        'uint8',   '8', ...
        'uint16',  '6', ...
        'uint32',  '2', ...
        'uint64',  '4' ...
    );
    a2acfg.itypes = { ...
        'int32', 'uint32', ...
        'int64', 'uint64', ...
        'int16', 'uint16', ...
        'int8',  'uint8',  ...
        'logical'};
    a2acfg.prcs = 8;
end

% number of arguments
if nargin < 1
    error('neuroelf:general:badNumberOfInputs', 'Invalid number of inputs.');
end

% setup vars
prcs = a2acfg.prcs;
dxct = 0;

% precision in second argument
if nargin > 1 && ...
    isnumeric(varargin{1}) && ...
   ~isempty(varargin{1})
    prcs = floor(varargin{1}(1));

% use hxdouble instead
elseif nargin > 1 && ...
    ischar(varargin{1}) && ...
    strcmpi(varargin{1}(:)', 'exact')
    dxct = 1;
end

% no negative precision values
if prcs < 0
    prcs = -prcs;
end

% if precision is 0, no decimal point
if prcs == 0
    prcstr = '%0.0f,';
    cpxstr = '%g + %gi,';

% otherwise use %g formatter
else
    prcstr = ['%.' int2str(prcs) 'g,'];
    cpxstr = [prcstr(1:end-1) ' + ' prcstr(1:end-1) 'i,'];
end

% non-builtin types would error, so those need overloading
try

    % try lookup from persistent car
    type = a2acfg.classes.(lower(class(anyvar)));

% handle errors
catch

    % error out if type unsupported
    error('neuroelf:any2ascii:invalidInputType', 'Invalid input type: %s.', class(anyvar));
end

% handle sparse matrices differently
if type == 'd' && ...
    issparse(anyvar)
    type = 'p';
end

% get dimensions
arraydims = size(anyvar);

% no array content
if any(arraydims == 0)

    % but not all dims zeros
    if ~all(arraydims == 0)

        % for numeric types
        if any('123456789cdip' == type)

            % simply build "empty" zeros with same size
            asciirep = ['[zeros(' any2ascii(arraydims) ')]'];

            % and prepend correct class cast
            if double(type) < 65
                asciirep = ['[' a2acfg.itypes{double(type) - 48} ...
                            '(' asciirep(2:(end-1))  ')]'];
            elseif type == 'c'
                asciirep = ['[char(' asciirep(2:(end-1)) ')]'];
            elseif type == 'i'
                asciirep = ['[single(' asciirep(2:(end-1)) ')]'];
            elseif type == 'p'
                asciirep = ['[sparse(' asciirep(2:(end-1)) ')]'];
            end

        % empty cell array
        elseif type == 'l'
            asciirep = ['[cell(' any2ascii(arraydims) ')]'];

        % empty struct
        else

            % get fieldnames
            fnams = fieldnames(anyvar);

            % empty field list
            if numel(fnams) == 0
                asciirep = ['[repmat(struct,' any2ascii(arraydims) ')]'];

            % or has fields
            else

                % create repmat/cellstruct call
                asciirep = sprintf('[repmat(cell2struct(cell([1,1,%d]),%s,3),%s)]', ...
                    numel(fnams), any2ascii(fnams'), any2ascii(arraydims));
            end
        end

    % all empty dims
    else

        % numeric type
        if any('123456789cdip' == type)

            % generate call
            if type == 'd'
                asciirep = '[]';
            elseif type == 'c'
                asciirep = '['''']';
            elseif type == 'i'
                asciirep = '[single([])]';
            elseif type == 'p'
                asciirep = '[sparse([])]';
            else
                asciirep = ['[' a2acfg.itypes{double(type) - 48} '([])]'];
            end

        % empty cell
        elseif type == 'l'
            asciirep = '{}';

        % empty struct
        else

            % get fieldnames
            fnams = fieldnames(anyvar);

            % empty field list
            if numel(fnams) == 0
                asciirep = ['[repmat(struct,' any2ascii(arraydims) ')]'];

            % or has fields
            else

                % create repmat/cellstruct call
                asciirep = sprintf('[repmat(cell2struct(cell([1,1,%d]),%s,3),%s)]', ...
                    numel(fnams), any2ascii(fnams'), any2ascii(arraydims));
            end
        end
    end
    return;
end

% full, 2D matrix
if type ~= 'p' && numel(arraydims) < 3

    % if this is not a cell array, use '[' brackets
    if type ~= 'l'
        asciirep = '[';

        % switch over type
        switch type

        % standard doubles/singles
        case {'d', 'i'}

            % doubles without hxdouble or complex content
            if type == 'd' && (~dxct || ~isreal(anyvar))

                % iterate over dims(1)
                numrep = cell(1, arraydims(1));
                for outer = 1:arraydims(1)

                    % try combined sprintf approach
                    if isreal(anyvar)
                        line = sprintf(prcstr, anyvar(outer, :));
                    else
                        line = sprintf(cpxstr, ...
                            lsqz([real(anyvar(outer, :)); imag(anyvar(outer, :))])');
                    end

                    % replace last comma by semicolon
                    line(end) = ';';

                    % put line into asciirep
                    numrep{outer} = line;
                end
                asciirep = [asciirep sprintf('%s', numrep{:})];

            % yet doubles (use hxdouble then)
            elseif type == 'd'

                % more "lines"
                if arraydims(1) ~= 1

                    % make a reshaped array
                    asciirep = ['[reshape(hxdouble(''' ...
                                  hxdouble(reshape(anyvar, 1, prod(arraydims))) ...
                                  '''),' ...
                                  sprintf('[%d,%d]', arraydims(1), arraydims(2)) ...
                                  ')]'];

                % for one-liners, it's OK to use hxdouble directly
                else
                    asciirep = ['[hxdouble(''' hxdouble(anyvar) ''')]'];
                end

            % only singles remain
            else

                % more "lines"
                if arraydims(1) ~= 1

                    % make a reshaped array
                    asciirep = ['[reshape(hxsingle(''' ...
                                  hxsingle(reshape(anyvar, 1, prod(arraydims))) ...
                                  '''),' ...
                                  sprintf('[%d,%d]', arraydims(1), arraydims(2)) ...
                                  ')]'];

                % for one-liners, it's OK to use hxsingle directly
                else
                    asciirep = ['[hxsingle(''' hxsingle(anyvar) ''')]'];
                end
            end

        % struct array, recursive call for each member value
        case {'s'}

            % get fieldnames and number of names
            fnames = fieldnames(anyvar);
            nnames = numel(fnames);

            % we have content
            if nnames > 0

                % just a single (1x1) struct
                if all(arraydims == 1)

                    % create struct call
                    asciirep = '[struct(';

                    % iterate over names
                    for fcount = 1:nnames

                        % get field (works only on 1x1)
                        cv = anyvar(1).(fnames{fcount});

                        % if not is cell, directly
                        if ~iscell(cv)
                            asciirep = [asciirep '''' fnames{fcount} ''',' ...
                                        any2ascii(cv, prcs, []) ',' ];

                        % otherwise within {} to get call correct
                        else
                            asciirep = [asciirep '''' fnames{fcount} ''',{' ...
                                        any2ascii(cv, prcs, []) '},' ];
                        end
                    end

                    % put to representation
                    asciirep = [asciirep(1:end-1) ');'];

                % multi-struct (array)
                else

                    asciirep = ['[cell2struct(' ...
                                  any2ascii(struct2cell(anyvar)) ',' ...
                                  any2ascii(fnames) ',1);'];
                end

            % no content (no fields, but size valid)
            else
                % just struct
                if all(arraydims == 1)
                    asciirep = [asciirep 'struct,'];

                % use cell2struct/struct2cell pair to convert
                else
                    asciirep = sprintf('[cell2struct(cell([0,%d,%d]),{},1)]', ...
                                         arraydims(1),arraydims(2));
                end
            end

        % character arrays
        case {'c'}

            % for empty character arrays
            if prod(arraydims) == 0

                % always give an empty array
                asciirep = [asciirep ''''','];

            % with content
            else

                % iterate over "lines"
                for outer = 1:arraydims(1)

                    % replace single quotes
                    ovar = strrep(anyvar(outer, :), '''', '''''');

                    % find unprintable characters
                    illegalc = find(ovar<32 | ovar>127);

                    % replace those
                    while ~isempty(illegalc)
                        ovar = strrep(ovar, ovar(illegalc(1)), ...
                               [''' char(' ...
                                sprintf('%d',double(ovar(illegalc(1)))) ...
                                ') ''']);

                        % find anew
                        illegalc = find(ovar<32 | ovar>127);

                        % but keep track of bracket usage!
                        prcs = -1;
                    end

                    % use brackets?
                    if prcs >= 0
                        asciirep = [asciirep '''' ovar '''' ';'];
                    else
                        asciirep = [asciirep '[''' ovar '''' '];'];
                    end
                end
            end

        % int/uint arrays
        otherwise

            % general case: everything but 1x1 logical
            if type ~= '9' || ...
                numel(anyvar) ~= 1
                asciirep = [asciirep class(anyvar) '(' ...
                            any2ascii(double(anyvar), varargin{2:end}) '),'];

            % 1x1 logical -> true
            elseif anyvar
                asciirep = [asciirep 'true,'];
            else
                asciirep = [asciirep 'false,'];
            end

        end

        % either add or replace closing brackets (upon content given)
        if asciirep(end) ~= '['
            asciirep(end) = ']';
        else
            asciirep(end + 1) = ']';
        end

        % third argument
        if nargin > 2 && ...
            isempty(varargin{2}) && ( ...
           ((type == 'd' || ...
             type == 's') && ...
             all(arraydims == 1)) || ...
            (type == 'c' && ...
             arraydims(1) == 1))

            % remove those brackets (for internal use)
            asciirep = asciirep(2:end-1);
        end

        % double brackets?
        while length(asciirep) > 2 && ...
            all(asciirep(1:2) == '[') && ...
            all(asciirep(end-1:end) == ']')

            % remove those anyway
            asciirep=asciirep(2:end-1);
        end

    % cell arrays
    else

        % correct opener
        asciirep = '{';

        % iterate over outer dimension
        for outer = 1:arraydims(1)

            % create "line"'s
            line = '';

            % iterate over inner dimension
            for inner = 1:arraydims(2)
                line = [line any2ascii(anyvar{outer,inner}, prcs, []) ','];
            end
            asciirep = [asciirep line(1:end-1) ';'];
        end

        % on non-empty content
        if asciirep(end) ~= '{'

            % replace final delimiter
            asciirep(end) = '}';

        % for empty content
        else

            % add closing bracket
            asciirep(end + 1) = '}';
        end
    end

% N > 2D, but not a sparse -> add reshape code
elseif type ~= 'p'

    % so pack the contents just as we suspect it to be packed :)
    origdims = sprintf('%d,', arraydims);
    asciirep = ['[reshape(' ...
                any2ascii(reshape(anyvar, 1, prod(arraydims)), varargin{2:end}) ...
                ',' origdims(1:(end-1)) ')]'];

% for sparse matrices, make special code
else

    % get content nicely done
    [spi, spj, sps] = find(anyvar);
    asciirep = ['[sparse(' any2ascii(spi, 0) ',' any2ascii(spj, 0) ',' ...
                any2ascii(sps,varargin{2:end}) ',' ...
                sprintf('%.0f,%.0f', arraydims(1), arraydims(2)) ')]'];

end

% lazy reformatted
if nargin > 2 && ...
  ~isempty(varargin{2}) && ...
   varargin{2} > 0
    lb = [' ...' char(10) char(9)];
    asciirep = strrep(strrep(strrep(asciirep, ...
                      ',',    ', '), ...
                      ';',   [';' lb]), ...
                      '), ', ['),' lb]);
end


% private function lsqz to linearly squeeze variable as in v(:)
function v = lsqz(v)
v = v(:);
