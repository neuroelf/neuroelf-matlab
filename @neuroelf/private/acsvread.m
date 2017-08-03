function [table, varargout] = acsvread(content, varargin)
% acsvread  - reads a text based (non-numeric) csv file
%
% due to the different approaches when reading csv files, both a cell
% and a struct oriented method can be used.
%
% FORMAT:         table = acsvread(file, [delimiter,] options)
%
% Input fields:
%    file         either filename (preferably absolute) to input or
%                 a newline seperated string containing a table or
%                 even a cell array with input lines to process
%    delimiter    delimiter character(s) (used for MATLAB's strtok)
%    options      struct array with multiple of the following options
%    .asmatrix    table array is a MxN cell array, not 1xM of 1xN
%    .convert     if given, try to convert numerical values to such;
%                 plus if set to 'deblank', deblanks all fields
%    .headline    if given, functions assumes named fields (heads);
%                 this results in a struct array, not a cell array!
%                 further, if headline is empty read headline from
%                 given input, otherwise use this at headline!
%    .multidelim  if given, treat delimiter as multiple tokens
%    .noquotes    if given, do not respect quotes (which is default)
%    .numonly     if given, only reads numbers
%    .range       if given, acsvread only reads within range can be
%                 either [upper-left-x,y, lower-right-x,y] coordinates
%                 or a string in MS-Excel notation ('A4:F9') where
%                 any missing coordinate is stretched to the edge
%    .readstart   if given, start reading at line after given token
%    .readstop    if given, stop reading at line before given token
%                 meaning: lines must contain this token to match
%    .replacef    cell array with "replace from" strings
%    .replacet    cell array with "replace to" string
%
% See also asciiread, splittocellc.

% Version:  v0.9b
% Build:    14012811
% Date:     Apr-08 2011, 10:18 PM EST
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

% enough arguments ?
if nargin < 1
    error( ...
        'neuroelf:TooFewArguments',...
        'Too few arguments. Try ''help %s''.',...
        mfilename ...
    );
end

% categorize arguments
table     = cell(0);
delimiter = ' ';
dodeblank = false;
opts      = struct;

if ischar(content)
    if any(content == 10)
        content(content == char(13)) = [];
    else
        content = asciiread(content, char(10));
    end
    content = splittocell(content, char(10), 0);
elseif ~iscell(content)
    error( ...
        'neuroelf:BadArgument',...
        'First argument must be a 1xN char or 1xN cell array!' ...
    );
end

if nargin > 1
    if nargin > 2
        if isstruct(varargin{2}), opts = varargin{2}; end
    end
    if isstruct(varargin{1})
        opts = varargin{1};
    else
        delimiter = sprintf(varargin{1});
    end
end
opts = opts(1);

if isfield(opts, 'asmatrix')
    asmatrix = 1;
else
    asmatrix = 0;
end
if isfield(opts, 'convert')
    convert = 1;
    if ischar(opts.convert)
        convertt = lower(opts.convert(:)');
        if numel(convertt) == 7 && ...
            all(convertt == 'deblank')
            dodeblank = true;
        end
    end
else
    convert = 0;
end
if isfield(opts, 'delimiter') && ...
    ischar(opts.delimiter)
    delimiter = opts.delimiter;
end
multidelim = isfield(opts, 'multidelim');
if isfield(opts, 'noquotes')
    quotes = false;
else
    quotes = true;
end
if isfield(opts, 'numonly')
    convert = 1;
    numonly = 1;
else
    numonly = 0;
end
if isfield(opts, 'range')
    range = opts.range;
    if ischar(range), range = exceltocoords(range); end
    if ~isnumeric(range)
        error( ...
            'neuroelf:BadOption',...
            'Range must be either numerical or in MS-Excel notation!' ...
        );
    end
else
    rs = '1';
    gt = 1;
    gm = size(content, 2);
    if isfield(opts, 'readstart')
        rt = opts.readstart;
        while gt <= gm && ...
            isempty(strfind(content{gt}, rt))
            gt = gt + 1;
        end
        gt = gt + 1;
        rs = num2str(gt);
    end
    rs = [rs ':'];
    if isfield(opts, 'readstop')
        rt = opts.readstop;
        while gt <= gm && ...
            isempty(strfind(content{gt}, rt))
            gt = gt + 1;
        end
        rs = [rs num2str(gt - 1)];
    else
        rs = [rs '65535'];
    end
    range = exceltocoords(rs);
end
if isfield(opts, 'replacef') && ...
    iscell(opts.replacef) && ...
    isfield(opts, 'replacet') && ...
    iscell(opts.replacet)
    replacef = opts.replacef;
    replacet = opts.replacet;
else
    replacef = cell(0);
    replacet = cell(0);
end

cs = size(content,2);
if range(2) > cs
    if nargout > 1
        varargout = cell(1, nargout - 1);
    end
    return;
end
if range(4) > cs || ...
   (range(4) == 65535 && ...
    cs > 65535)
    range(4) = cs;
end
content = content(range(2):range(4));
rcs     = size(content, 2);

if numel(replacef) > 0 && ...
    numel(replacet) == numel(replacef)
    for rn = 1:numel(replacef)
        for lc = 1:rcs
            content{lc} = strrep(content{lc}, replacef{rn}, replacet{rn});
        end
    end
end

headline = '';
if isfield(opts, 'headline')
    headline = opts.headline;
    if ~ischar(headline) || ...
        isempty(headline)
        if rcs
            headline = content{1};
            content  = content(2:rcs);
            rcs      = rcs - 1;
        else
            error( ...
                'neuroelf:acsvread:EmptyTable',...
                'Empty table insuitable for empty headline.' ...
            );
        end
    else
        headline = sprintf(headline);
    end
end
cs = size(content, 2);

%  - work on headline if needed and choose for what to do...
if size(headline, 2) > 0
    tokens    = splittocellc(headline, delimiter, multidelim, multidelim);
    numtokens = size(tokens, 2);
    for ntdb  = 1:numtokens
        tokens{ntdb} = deblank(tokens{ntdb});
    end
else
    numtokens = 0;
end

if numtokens == 0
    table   = cell(1,cs);
else
    xtokens = cell(0);
    ttokens = struct;
    for tc = 1:numtokens
        ttoken = makelabel(tokens{tc});
        if isfield(ttokens, ttoken)
            ttoken = sprintf('%s_%d', ttoken, tc);
        end
        ttokens.(ttoken) = 1;
        dc = 2 * tc;
        xtokens((dc - 1):dc) = {ttoken, ''};
    end
    table = struct(xtokens{:});
    table(cs) = table(1);
end
if rcs == 0
    if nargout > 1
        varargout{1} = 0;
        if nargout > 2
            varargout{2} = 0;
        end
    end
    return;
end

rfrom = range(1);
rto   = range(3);
if numtokens == 0
    rmaxsize = 0;
    for lc = 1:cs
        if convert == 0
            table{lc} = ...
                splittocellc(content{lc}, delimiter, multidelim, multidelim, quotes);
        else
            table{lc} = i_convertposs( ...
                splittocellc(content{lc}, delimiter, multidelim, multidelim, quotes), ...
                dodeblank);
        end
        rmaxsize = max(rmaxsize, size(table{lc},2));
    end
    if rfrom > rmaxsize
        error( ...
            'neuroelf:OutOfRangeError',...
            'Range for table is invalid.' ...
        );
    end
    if rto   > rmaxsize, rto = rmaxsize; end
    for lc = 1:cs
        rrealsize = size(table{lc}, 2);
        if rrealsize < rmaxsize
            [table{lc}{(rrealsize + 1):rmaxsize}] = deal('');
        end
        if rfrom > 1 || ...
            rto < rmaxsize
            table{lc} = table{lc}(rfrom:rto);
        end
    end
    if asmatrix
        mtable = cell(length(table), length(table{1}));
        for lc = 1:cs
            mtable(lc, :) = table{lc};
        end
        table = mtable;
    end
else
    ntt = 2:2:(2 * numtokens);
    for lc = 1:cs
        if convert == 0
            values = ...
                splittocellc(content{lc}, delimiter, multidelim, multidelim, quotes);
        else
            values = i_convertposs( ...
                splittocellc(content{lc}, delimiter, multidelim, multidelim, quotes), ...
                dodeblank);
        end
        numvalues = size(values, 2);
        if numvalues < numtokens
            [values{(numvalues+1):numtokens}] = deal('');
        end
        xtokens(ntt) = values(1:numel(ntt));
        try
            table(lc) = struct(xtokens{:});
        catch ne_eo;
            warning( ...
                'neuroelf:ConversionFailed', ...
                'Couldn''t convert CSV line %d:%s%s%s%s', ...
                lc, char(10), ...
                content{lc}(1:min(length(content{lc}), 80)), ...
                char(10), ne_eo.message ...
            );
        end
    end
end
if numonly
    tlinen = true(1, numel(table{1}));
    for lc = 1:numel(table{1})
        tlinen(lc) = isnumeric(table{1}{lc}) && numel(table{1}{lc}) == 1;
    end
    tlinen = find(tlinen);
    ntable = zeros(numel(table), numel(tlinen));
    for lc = 1:numel(table)
        for cc = 1:numel(tlinen)
            try
                ntable(lc, cc) = table{lc}{tlinen(cc)};
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
                ntable(lc, cc) = NaN;
            end
        end
    end
    table = ntable;
end

if nargout > 1
    if numtokens == 0
        varargout{1} = size(table, 1);
    else
        varargout{1} = size(table, 2);
    end
    if nargout > 2
        if numtokens == 0
            varargout{2} = rmaxsize;
        else
            varargout{2} = numtokens;
        end
    end
end

return;


%%% internal functions


function values = i_convertposs(values, dodeblank)
for vc = 1:size(values, 2)
    if numel(values{vc}) > 0
        try
            valdbl = u8str2double(values{vc}, 1, 1);
            if ~isnan(valdbl) && ...
               ~isinf(valdbl)
                values{vc} = valdbl;
                continue;
            end
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
        end
        valdbl = deblank(values{vc});
        if ~isempty(valdbl) && ...
           ~dodeblank
            if valdbl(1) == '$'
                try
                    values{vc} = eval(valdbl(2:end));
                catch ne_eo;
                    neuroelf_lasterr(ne_eo);
                    disp(['acsvread> Error parsing numeric value: ' valdbl(2:end)]);
                end
            elseif values{vc}(1) == ']'
                try
                    values{vc} = evalin('base', valdbl(2:end));
                catch ne_eo;
                    neuroelf_lasterr(ne_eo);
                    disp(['acsvread> Error evaluating expression: ' valdbl(2:end)]);
                end
            elseif strcmpi(valdbl, 'inf')
                values{vc} = Inf;
            elseif strcmpi(valdbl, '-inf')
                values{vc} = -Inf;
            elseif strcmpi(valdbl, 'nan')
                values{vc} = NaN;
            else
                values{vc} = valdbl;
            end
        elseif dodeblank || isempty(valdbl)
            values{vc} = valdbl;
        end
    end
end
