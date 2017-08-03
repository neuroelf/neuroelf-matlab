function plp = importplp(matfile, opts)
% importplp  - import a PLP object (from a MKDA DB)
%
% FORMAT:       plp = importplp(DBfile [, opts])
%
% Input fields:
%
%       DBfile      MAT file containing the MKDA DB
%       opts        optional settings
%        .discard   cell array with columns to discard (default: {})
%        .icase     ignore-case flag for PLP::AddColumn call (true)
%        .regexprep Sx2 cell array with regexprep replacements (from, to)
%        .strrep    Sx2 cell array with strrep replacements (from, to)
%        .studycol  name of study column (defaut: 'study')
%        .tal2x     either of {'icbm'}, 'mni', or 'none'
%        .talcode   token recognized as Talairach, default 'T88'
%        .valcode   Vx2 cell array with value codings ({'token', [value]}
%                   the default is to code no = 0, yes = 1
%
% Output fields:
%
%       plp         PLP object
%
% Note: the replacements are done in the order 1.) regexprep, 2.) strrep
%
% Note: any text-based coordinate database MUST contain the following
%       columns (case insensitive): Study, X, Y, Z; in addition, the name
%       of a coordinate system convention column, if present, must begin
%       with "coor" and must NOT contain "orig" and all cells/values in
%       this column must be of type char and non-empty.

% Version:  v1.1
% Build:    16051821
% Date:     May-18 2016, 9:00 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011 - 2015, 2016, Jochen Weber
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

% argument check
if nargin < 1 || ~ischar(matfile) || isempty(matfile)
    error('neuroelf:general:badArgument', 'Bad or missing argument.');
end
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'discard') || ~iscell(opts.discard) || isempty(opts.discard)
    opts.discard = {};
else
    opts.discard = opts.discard(:);
    for ec = numel(opts.discard):-1:1
        if ~ischar(opts.discard{ec}) || isempty(opts.discard{ec}) || ...
            numel(opts.discard{ec}) ~= size(opts.discard{ec}, 2) || ...
           ~isvarname(opts.discard{ec})
            opts.discard(ec) = [];
        end
    end
    opts.discard = unique(opts.discard)';
end
if ~isfield(opts, 'icase') || ~islogical(opts.icase) || numel(opts.icase) ~= 1
    opts.icase = true;
end
if ~isfield(opts, 'regexprep') || ~iscell(opts.regexprep) || ...
    ndims(opts.regexprep) > 2 || size(opts.regexprep, 2) ~= 2
    opts.regexprep = cell(0, 2);
else
    for ec = size(opts.regexprep, 1):-1:1
        if ~ischar(opts.regexprep{ec, 1}) || isempty(opts.regexprep{ec, 1}) || ...
            numel(opts.regexprep{ec, 1}) ~= size(opts.regexprep{ec, 1}, 2) || ...
           ~ischar(opts.regexprep{ec, 2})
            opts.regexprep(ec, :) = [];
        else
            if isempty(opts.regexprep{ec, 2})
                opts.regexprep{ec, 2} = '';
            else
                opts.regexprep{ec, 2} = opts.regexprep{ec, 2}(:)';
            end
        end
    end
end
if ~isfield(opts, 'strrep') || ~iscell(opts.strrep) || ...
    ndims(opts.strrep) > 2 || size(opts.strrep, 2) ~= 2
    opts.strrep = cell(0, 2);
else
    for ec = size(opts.strrep, 1):-1:1
        if ~ischar(opts.strrep{ec, 1}) || isempty(opts.strrep{ec, 1}) || ...
            numel(opts.strrep{ec, 1}) ~= size(opts.strrep{ec, 1}, 2) || ...
           ~ischar(opts.strrep{ec, 2})
            opts.strrep(ec, :) = [];
        else
            if isempty(opts.strrep{ec, 2})
                opts.strrep{ec, 2} = '';
            else
                opts.strrep{ec, 2} = opts.strrep{ec, 2}(:)';
            end
        end
    end
end
if ~isfield(opts, 'studycol') || ~ischar(opts.studycol) || ...
    isempty(opts.studycol) || ~isrealvarname(opts.studycol(:)')
    opts.studycol = 'study';
end
if ~isfield(opts, 'tal2x') || ~ischar(opts.tal2x) || ...
    isempty(opts.tal2x) || ~any(lower(opts.tal2x(1)) == 'imn')
    if isempty(regexpi(matfile(:)', '\.mat$'))
        opts.tal2x = 'i';
    else
        opts.tal2x = 'n';
    end
else
    opts.tal2x = lower(opts.tal2x(1));
end
if ~isfield(opts, 'talcode') || ~ischar(opts.talcode) || isempty(opts.talcode)
    opts.talcode = 'T88';
else
    opts.talcode = opts.talcode(:)';
end
if ~isfield(opts, 'valcode') || ~iscell(opts.valcode) || ...
    ndims(opts.valcode) > 2 || size(opts.valcode, 2) ~= 2
    opts.valcode = {'no', 0; 'yes', 1};
else
    for ec = size(opts.valcode, 1):-1:1
        if ~ischar(opts.valcode{ec, 1}) || isempty(opts.valcode{ec, 1}) || ...
            numel(opts.valcode{ec, 1}) ~= size(opts.valcode{ec, 1}, 2) || ...
           ~isa(opts.valcode{ec, 2}, 'double') || numel(opts.valcode{ec, 2}) ~= 1
            opts.valcode(ec, :) = [];
        end
    end
end

% try to load DB
matfile = matfile(:)';
if ~any(matfile == '/' | matfile == '\')
    matfile = [pwd '/' matfile];
end

% some error might well happen...
try
    % preset to detect xls/txt files
    dbh = [];

    % for EXCEL files
    if ~isempty(regexpi(matfile, '\.xlsx?$'))

        % read in content
        [dbh, db, dbc] = xlsread(matfile);

        % remove empty header columns
        for cc = size(dbc, 2):-1:1
            if ~isempty(dbc{1, cc})
                break;
            end
            dbc(:, cc) = [];
        end

        % set header
        dbh = dbc(1, :);

    % for text files
    elseif ~isempty(regexpi(matfile, '\.(csv|txt)$'))
        db = splittocellc(asciiread(matfile), char([10, 13]), true, true);
        if any(db{1} == char(9))
            fdelim = char(9);
        elseif any(db{1} == '|')
            fdelim = '|';
        elseif any(db{1} == ';')
            fdelim = ';';
        elseif any(db{1} == ',')
            fdelim = ',';
        else
            error( ...
                'neuroelf:UnsupportedFormat', ...
                'Text file not supported.' ...
            );
        end

        % read by default without heeding quotes
        dbc = acsvread(matfile, fdelim, ...
            struct('asmatrix', 1, 'convert', 1, 'noquotes', 1));

        % test if last column has too many missing values
        if sum(cellfun('isempty', dbc(:, end))) > (0.9 * size(dbc, 1))

            % re-read with taking quotes into account
            dbc = acsvread(matfile, fdelim, struct('asmatrix', 1, 'convert', 1));
        end

        % check content
        if isa(dbc{1}, 'double') && numel(dbc{1}) == 1 && ...
            abs(dbc{1} - size(dbc, 2)) <= 2
            dbc = dbc(2:end, :);
        end
        dbh = dbc(1, :)';
        for cc = numel(dbh):-1:1
            if ~isempty(dbh{cc})
                break;
            end
            dbh(cc) = [];
        end
    end

    % content read already
    if ~isempty(dbh)
        if ~any(strcmpi(opts.studycol, dbh)) || ~any(strcmpi('x', dbh)) || ...
           ~any(strcmpi('y', dbh)) || ~any(strcmpi('z', dbh)) || ...
            any(strcmpi('', dbh))
            error('neuroelf:import:unsupportedFormat', ...
                'Crucial column missing or empty column name.');
        end
        for cc = size(dbc, 1):-1:1
            if all(cellfun('isempty', dbc(cc, :)))
                dbc(cc, :) = [];
            end
        end

        % create database
        db = struct('DB', struct);
        for fc = 1:numel(dbh)
            dbh{fc} = makelabel(dbh{fc});
            if numel(dbh{fc}) == 1 && ...
                any('XYZ' == dbh{fc})
                dbh{fc} = lower(dbh{fc});
            elseif strcmpi(dbh{fc}, opts.studycol)
                dbh{fc} = 'Study';
            end
            if isfield(db.DB, dbh{fc})
                for sc = 2:64
                    if ~isfield(db.DB, sprintf('%s%02d', dbh{fc}, sc))
                        dbh{fc} = sprintf('%s%02d', dbh{fc}, sc);
                        break;
                    end
                end
                if isfield(db.DB, dbh{fc})
                    error('neuroelf:import:inputError', ...
                        'Too many columns with the same name.');
                end
            end
            if ~any(cellfun(@ischar, dbc(2:end, fc)))
                db.DB.(dbh{fc}) = cat(1, dbc{2:end, fc});
            else
                db.DB.(dbh{fc}) = dbc(2:end, fc);
            end
        end

    % otherwise, try as MAT file
    else
        db = load(matfile, '-mat');
    end

    % now check the format
    if ~isfield(db, 'DB') || ~isstruct(db.DB) || numel(db.DB) ~= 1
        error('neuroelf:import:badArgument', 'Bad structure of DB matfile.');
    end
    if ~isfield(db.DB, 'Study') || ~isfield(db.DB, 'x') || ...
       ~isfield(db.DB, 'y') || ~isfield(db.DB, 'z')
        error('neuroelf:import:badArgument', 'Missing fields in DB.');
    end
    if iscell(db.DB.Study) && any(cellfun('isempty', db.DB.Study))
        error('neuroelf:import:badArgument', ...
            'Study field must not be empty for any data in DB.');
    end
    if ~isa(db.DB.x, 'double') || numel(db.DB.x) ~= size(db.DB.x, 1) || ...
       ~isa(db.DB.y, 'double') || ~isequal(size(db.DB.y), size(db.DB.x)) || ...
       ~isa(db.DB.z, 'double') || ~isequal(size(db.DB.z), size(db.DB.x))
        error('neuroelf:import:badArgument', ...
            'Invalid coordinate data in DB (please check for additional commas, etc.).');
    end
catch ne_eo;
    rethrow(ne_eo);
end
db = db.DB;

% get fieldnames
dbf = fieldnames(db);

% provide support for T88 coordinate system ala TW tools
if opts.tal2x ~= 'n'
    csysc = find(~cellfun('isempty', regexpi(dbf, '^coor')));
    if isempty(csysc)
        csysc = find(~cellfun('isempty', regexpi(dbf, '^space$')));
    end
else
    csysc = [];
end
if ~isempty(csysc)
    dbc = struct2cell(db);
    xcol = findfirst(strcmpi(dbf, 'x'));
    ycol = findfirst(strcmpi(dbf, 'y'));
    zcol = findfirst(strcmpi(dbf, 'z'));
end
for cc = 1:numel(csysc)
    if isempty(regexpi(dbf{csysc(cc)}, 'ori?g(inal)?$'))
        if ~all(cellfun(@ischar, dbc{csysc(cc)})) || ...
           ~any(strcmpi(dbc{csysc(cc)}, opts.talcode))
            continue;
        end
        dbnc = dbc([csysc(cc), xcol, ycol, zcol]);
        dbcc = false(numel(dbnc{1}), 1);
        for rc = 2:numel(dbcc)
            if strcmpi(dbnc{1}{rc}, opts.talcode)
                dbcc(rc) = true;
                try
                    xyz = cat(2, dbnc{2}(rc), dbnc{3}(rc), dbnc{4}(rc));
                    if ~isequal(size(xyz), [1, 3])
                        error('Bad coordinate content!');
                    end
                    if opts.tal2x == 'i'
                        xyz = tal2icbm(xyz);
                    else
                        xyz = tal2mni(xyz);
                    end
                    dbnc{1}{rc} = 'icbm';
                    dbnc{2}(rc) = xyz(1);
                    dbnc{3}(rc) = xyz(2);
                    dbnc{4}(rc) = xyz(3);
                catch ne_eo;
                    neuroelf_lasterr(ne_eo);
                    warning( ...
                        'neuroelf:ConversionError', ...
                        'Error converting coordinate %d: %s.', ...
                        rc - 1, ne_eo.message ...
                    );
                    dbcc(:) = false;
                    break;
                end
            end
        end
        db.(dbf{csysc(cc)}) = dbnc{1};
        db.(dbf{xcol}) = dbnc{2};
        db.(dbf{ycol}) = dbnc{3};
        db.(dbf{zcol}) = dbnc{4};
        if opts.tal2x == 'i'
            opts.tal2x = 'ICBM';
        else
            opts.tal2x = 'MNI';
        end
        disp(sprintf( ...
            'Converted %d coordinates from %s to %s.', sum(dbcc), ...
            opts.talcode, opts.tal2x));
        break;
    end
end

% replace contents
dbft = emptystruct(dbf(:), [1, 1]);
for fc = 1:numel(dbf)
    dbft.(dbf{fc}) = false;
    if iscell(db.(dbf{fc}))
        if size(db.(dbf{fc}), 2) ~= 1 || ...
            ndims(db.(dbf{fc})) > 2
            error( ...
                'neuroelf:BadArgument', ...
                'Invalid cell array content for field %s in DB.', ...
                dbf{fc} ...
            );
        end
        charcell = cellfun(@ischar, db.(dbf{fc}));
        if sum(charcell) >= (0.5 * numel(db.(dbf{fc})))
            dbft.(dbf{fc}) = true;
        end
        try
            for rc = 1:size(opts.regexprep, 1)
                db.(dbf{fc})(charcell) = regexprep(db.(dbf{fc})(charcell), ...
                    opts.regexprep{rc, 1}, opts.regexprep{rc, 2});
            end
        catch ne_eo;
            error( ...
                'neuroelf:BadArgument', ...
                'Error applying regexprep (''%s'' -> ''%s''): %s.', ...
                opts.regexprep{rc, 1}, opts.regexprep{rc, 2}, ne_eo.message ...
            );
        end
        try
            for rc = 1:size(opts.strrep, 1)
                db.(dbf{fc})(charcell) = strrep(db.(dbf{fc})(charcell), ...
                    opts.strrep{rc, 1}, opts.strrep{rc, 2});
            end
        catch ne_eo;
            error( ...
                'neuroelf:BadArgument', ...
                'Error applying strrep (''%s'' -> ''%s''): %s.', ...
                opts.strrep{rc, 1}, opts.strrep{rc, 2}, ne_eo.message ...
            );
        end
        try
            for rc = 1:size(opts.valcode, 1)
                charcell = find(cellfun(@ischar, db.(dbf{fc})));
                charcell(~strcmpi(db.(dbf{fc})(charcell), opts.valcode{rc, 1})) = [];
                if ~isempty(charcell)
                    db.(dbf{fc})(charcell) = ...
                        repmat(opts.valcode(rc, 2), numel(charcell), 1);
                end
            end
            charcell = cellfun(@ischar, db.(dbf{fc}));
            if ~any(charcell)
                db.(dbf{fc}) = cat(1, db.(dbf{fc}));
                if numel(db.(dbf{fc})) ~= numel(charcell)
                    error( ...
                        'neuroelf:BadArgument', ...
                        'Badly formated cells in field %s', ...
                        dbf{fc} ...
                    );
                end
            end
        catch ne_eo;
            error( ...
                'neuroelf:BadArgument', ...
                'Error replacing ''%s'' with value %g: %s.', ...
                opts.valcode{rc, 1}, opts.valcode{rc, 2}, ne_eo.message ...
            );
        end
        charcell = cellfun(@ischar, db.(dbf{fc}));
        if sum(charcell) < (0.5 * numel(db.(dbf{fc})))
            dbft.(dbf{fc}) = false;
        end
    end
end

% remove unwanted
opts.discard = [opts.discard; {'x'; 'y'; 'z'; 'Study'}];
unw = multimatch(dbf, opts.discard);
dbf(unw > 0) = [];

% create plp
plp = xff('new:plp');

% add default columns at the beginning
plp.AddColumn('X', db.x);
plp.AddColumn('Y', db.y);
plp.AddColumn('Z', db.z);
plp.AddColumn('Study', db.Study, opts.icase);

% then add all remaining columns
plpct = struct('X', false, 'Y', false, 'Z', false, 'Study', dbft.Study);
for fc = 1:numel(dbf)
    plp.AddColumn(dbf{fc}, db.(dbf{fc}), opts.icase);
    plpct.(dbf{fc}) = dbft.(dbf{fc});
end

% by default, use all points
upc = findfirst(strcmpi(plp.ColumnNames, 'usepoint'));
if ~isempty(upc) && ~any(strcmpi(fieldnames(db), 'usepoint'))
    plp.Points(:, upc) = 1;
end

% remove double quotes from labels
l = plp.Labels;
for lc = 1:numel(l)
    l{lc}(l{lc} == '"') = [];
end
plp.Labels = l;

% remove empty columns
ecols = find(all(plp.Points == 0, 1));
ecoln = plp.ColumnNames(ecols);
plp.ColumnNames(ecols) = [];
plp.Points(:, ecols) = [];
for ecolc = 1:numel(ecoln)
    if isfield(plpct, ecoln{ecolc})
        plpct = rmfield(plpct, ecoln{ecolc});
    end
end
plp.NrOfColumns = numel(plp.ColumnNames);

% add RunTimeVars
plp.RunTimeVars.AutoSave = true;
plp.RunTimeVars.SourceFile = matfile;
plp.RunTimeVars.ColumnIsText = plpct;
