function epti = readeprimetextlog(filename, flat)
% readeprimetextlog  - read eprime text log file
%
% FORMAT:       epti = readeprimetextlog(filename [, flat])
%
% Input fields:
%
%       filename    filename of text-based logfile
%       flat        if given and true, reduce all levels to single level
%
% Output fields:
%
%       epti        eprime text info

% Version:  v1.0
% Build:    16022509
% Date:     Feb-25 2016, 9:32 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011 - 2016, Jochen Weber
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
if nargin < 1 || isempty(filename) || (~iscell(filename) && ...
    (~ischar(filename) || exist(regexprep(filename(:)', '\,\d+$', ''), 'file') ~= 2))
    error('neuroelf:BadArgument', 'Bad or missing argument.');
end
if ischar(filename)
    filename = filename(:)';
    filecont = '';
else
    filecont = filename;
    filename = 'cellarray';
end

% read content
try
    if ~isempty(regexpi(filename, '\.xlsx?(\,\d+)?$'))
        hassheet = findfirst(filename == ',', -1);
        if isempty(hassheet)
            [headers, table, filecont] = xlsread(filename);
        else
            [headers, table, filecont] = xlsread(filename(1:hassheet-1), ...
                str2double(filename(hassheet+1:end)));
        end
    end
    if ~isempty(filecont)
        filecont(cellfun('isempty', filecont)) = {''};
        filedouble = cellfun('isclass', filecont, 'double');
        for lc = 1:numel(filecont)
            if filedouble(lc)
                filecont{lc} = sprintf('%g', filecont{lc});
            end
        end
        for lc = 1:size(filecont, 1)
            filecont{lc, 1} = gluetostring(filecont(lc, :), char(9));
        end
        filecont = filecont(:, 1);
    end
catch ne_eo;
    neuroelf_lasterr(ne_eo);
end
if isempty(filecont)
    try
        filecont = splittocellc(asciiread(filename), char([10, 13]), true, true);
        filecont = filecont(:);
    catch ne_eo;
        rethrow(ne_eo);
    end
end

% tabular file guess
if numel(filecont) > 2 && ...
    any(~cellfun('isempty', regexpi(filecont(1:3), '^experiment')))

    % skip first line if with filename, and then any empty lines
    if ~isempty(strfind(filecont{1}, '.edat'))
        filecont(1) = [];
    end
    filecont(cellfun('isempty', filecont)) = [];

    % several "Experiment" lines
    elines = find(~cellfun('isempty', regexpi(filecont, '^experiment')));
    if numel(elines) > 1
        try
            epti = readeprimetextlog([{'.edat'}; filecont(1:(elines(2)-1))]);
            elines(end+1) = numel(filecont) + 1;
            for cc = 2:(numel(elines)-1)
                pepti = readeprimetextlog([{'.edat'}; ...
                    filecont(elines(cc):(elines(cc+1)-1))]);
                epti.Log = catstruct(epti.Log, pepti.Log);
            end
        catch ne_eo;
            rethrow(ne_eo);
        end
        return;
    end

    % parse rest into cell arrays
    headers = splittocellc(filecont{1}, char(9));
    headers = headers(:);
    table = cell(numel(headers), numel(filecont) - 1);
    for lc = 2:numel(filecont)
        rowcont = splittocellc(filecont{lc}, char(9));
        if numel(rowcont) < numel(headers)
            table(1:numel(rowcont), lc - 1) = rowcont(:);
        else
            table(1:numel(headers), lc - 1) = lsqueeze(rowcont(1:numel(headers)));
        end
    end

    % replace numbers
    tnumber = find(~cellfun('isempty', regexpi(table(:), ...
        '^[\+\-]?(\d+|\d+\.\d+|\.\d+)([eE][\+\-]?\d+)?$')));
    for lc = 1:numel(tnumber)
        table{tnumber(lc)} = str2double(table{tnumber(lc)});
    end

    % make sure headers are labels
    for lc = 1:numel(headers)
        headers{lc} = makelabel(headers{lc});
    end
    if numel(unique(headers)) ~= numel(headers)
        error( ...
            'neuroelf:BadFileContent', ...
            'Header fields must be unique.' ...
        );
    end

    % create struct
    tstruct = cell2struct(table, headers, 1);

    % create output
    epti = struct( ...
        'Experiment',          'unknown', ...
        'Display_RefreshRate', 60, ...
        'Group',               1, ...
        'LevelName',           {{'FlatTable'}}, ...
        'Log',                 {{struct}}, ...
        'RandomSeed',          round((2^32 - 1) * (rand(1, 1) - 0.5)), ...
        'SessionDate',         datestr(now, 'mm-dd-yyyy'), ...
        'SessionTime',         datestr(now, 13), ...
        'Session',             1, ...
        'Subject',             1, ...
        'VersionPersist',      1);

    % extract *some* fields
    if isfield(tstruct, 'Experiment')
        epti.Experiment = tstruct(1).Experiment;
    elseif isfield(tstruct, 'ExperimentName')
        epti.Experiment = tstruct(1).ExperimentName;
    end
    if isfield(tstruct, 'Display_RefreshRate')
        epti.Display_RefreshRate = tstruct(1).Display_RefreshRate;
    end
    if isfield(tstruct, 'Group')
        epti.Group = tstruct(1).Group;
    end
    if isfield(tstruct, 'RandomSeed')
        epti.RandomSeed = tstruct(1).RandomSeed;
    end
    if isfield(tstruct, 'Session')
        epti.Session = tstruct(1).Session;
    end
    if isfield(tstruct, 'SessionDate')
        epti.SessionDate = tstruct(1).SessionDate;
    end
    if isfield(tstruct, 'SessionTime')
        epti.SessionTime = tstruct(1).SessionTime;
    end
    if isfield(tstruct, 'Subject')
        epti.Subject = tstruct(1).Subject;
    end

    % set in epti
    epti.Log = tstruct;

    % return
    return;
end

% remove empty lines
for lc = 1:numel(filecont)
    filecont{lc} = ddeblank(filecont{lc});
end
filecont(cellfun('isempty', filecont)) = [];

% does not contains any "Level" tags
lg = grep(filecont, 'Level');
if isempty(lg)
    
    % parsing as simple table
    if any(filecont{1} == char(9))
        fsep = {char(9), false, true};
    elseif any(filecont{1} == ',')
        fsep = {',', false, true};
    else
        fsep = {' ', true, true};
    end
    logfields = splittocellc(filecont{1}, fsep{:});
    while ~isempty(logfields) && ...
        isempty(logfields{end})
        logfields(end) = [];
    end
    for lc = 1:numel(logfields)
        logfields{lc} = makelabel(logfields{lc});
        if lc > 1 && ...
            any(strcmpi(logfields(1:lc-1), logfields{lc}))
            logfields{lc} = sprintf('%s_%d', logfields{lc}, lc);
        end
    end
    if sum(cellfun('isempty', regexp(logfields, '^V_'))) > (0.5 * numel(logfields))
        filecont(1) = [];
    else
        for lc = 1:numel(logfields)
            logfields{lc} = sprintf('column%03d', lc);
        end
    end
    filecont = filecont(:);
    filecont{1, numel(logfields)} = filecont{1};
    for lc = 1:size(filecont, 1)
        tablefields = splittocellc(filecont{lc, 1}, fsep{:});
        tablefields = tablefields(:)';
        ntf = min(numel(logfields), numel(tablefields));
        filecont(lc, 1:ntf) = tablefields(1, 1:ntf);
    end
    filecont(cellfun('isempty', filecont)) = {''};
    filenum = ~cellfun('isempty', regexpi(filecont, '^\d+(\.\d+)?([eE][\+\-]\d+)?$')) | ...
       ~cellfun('isempty', regexpi(filecont, '^[\+\-]?(inf|nan)$'));
    for lc = 1:numel(filecont)
        if filenum(lc)
            filecont{lc} = str2double(filecont{lc});
        end
    end
    
    % create output
    epti = struct( ...
        'Experiment',          'unknown', ...
        'Display_RefreshRate', 60, ...
        'Group',               1, ...
        'LevelName',           {{'FlatTable'}}, ...
        'Log',                 cell2struct(filecont, logfields(:), 2), ...
        'RandomSeed',          round((2^32 - 1) * (rand(1, 1) - 0.5)), ...
        'SessionDate',         datestr(now, 'mm-dd-yyyy'), ...
        'SessionTime',         datestr(now, 13), ...
        'Session',             1, ...
        'Subject',             1, ...
        'VersionPersist',      1);
    return;
end

% contains header
if ~isempty(regexpi(filecont{1}, 'header\s+start'))

    % split off header
    for lc = 2:numel(filecont)
        if ~isempty(regexpi(filecont{lc}, 'header\s+end'))
            break;
        end
    end
    epti = rep_parseframe([{'Level: 1'}; filecont(1:lc)]);
    epti.Log = {struct};
    filecont = filecont((lc+1):end);
else
    epti = struct( ...
        'Experiment',          'unknown', ...
        'Display_RefreshRate', 60, ...
        'Group',               1, ...
        'LevelName',           {{'Session', 'Block', 'Trial', 'SubTrial', 'LogLevel5'}}, ...
        'Log',                 {{struct}}, ...
        'RandomSeed',          round((2^32 - 1) * (rand(1, 1) - 0.5)), ...
        'SessionDate',         datestr(now, 'mm-dd-yyyy'), ...
        'SessionTime',         datestr(now, 13), ...
        'Session',             1, ...
        'Subject',             1, ...
        'VersionPersist',      1);
end

% find level lines
ll = find(~cellfun('isempty', regexpi(filecont(1:end-1), '^level\:\s+\d+$')) & ...
    ~cellfun('isempty', regexpi(filecont(2:end), 'logframe\s+start')));
ll(end+1) = numel(filecont) + 1;

% parse frames
f = cell(1, 10);
fl = 0;
for lc = 1:(numel(ll) - 1);

    % parse frame and get level
    try
        [fr, l] = rep_parseframe(filecont(ll(lc):(ll(lc+1)-1)));
    catch ne_eo;
        rethrow(ne_eo);
    end

    % stack previous frames
    if l < fl
        if l ~= (fl - 1)
            error( ...
                'neuroelf:BadHierarchy', ...
                'Bad log frame hierarchy detected.' ...
            );
        end
        fr.(sprintf('Level%d', fl)) = f{fl};
        f{fl} = [];
    end
    fl = l;

    % put into correct level
    if isempty(f{fl})
        f{fl} = fr;
    else
        f{fl} = catstruct(f{fl}(:), fr);
    end
end

% test last level
if fl ~= 1
    error( ...
        'neuroelf:ProcessError', ...
        'Error processing file: last level must be 1.' ...
    );
end

% stack into header
epti.Log = f{1};

% flat
if nargin > 1 && ...
    islogical(flat) && ...
    numel(flat) == 1 && ...
    flat
    epti.Log = flatlog(epti.Log, 2);
end



% sub-function: rep_parseframe
function [f, lv] = rep_parseframe(l)

% too few lines or bad content
if numel(l) < 3 || ...
    isempty(regexpi(l{2}, '\s+start\s+')) || ...
    isempty(regexpi(l{end}, '\s+end\s+'))
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid log frame layout.' ...
    );
end

% first line must contain level
if isempty(regexpi(l{1}, '^level\:\s+\d+$'))
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid level indication in frame header: %s.', ...
        l{1} ...
    );
end
lv = regexprep(l{1}, 'level\:\s+', '', 'ignorecase');
lv = str2double(lv);

% parse the rest of the lines into a struct
f = struct;
for lc = 3:(numel(l) - 1)

    % get fieldname and value
    [fn, fv] = strtok(l{lc}, ':');
    fn = makelabel(fn);
    fv = ddeblank(fv(2:end));

    % value numerical
    if ~isempty(regexpi(fv, '^[\+\-]?(\d+|\d+\.\d+|\.\d+)([eE][\+\-]?\d+)?$'))
        fv = str2double(fv);
    end

    % value doesn't exist
    if ~isfield(f, fn)

        % assign
        f.(fn) = fv;

    % otherwise
    else

        % all numbers
        if isa(fv, 'double') && ...
            isa(f.(fn), 'double')

            % put at the end
            f.(fn) = [f.(fn)(:)', fv];

        % or
        else

            % convert to cell if necessary
            if ~iscell(f.(fn))
                f.(fn) = {f.(fn)};
            end

            % add to the end
            f.(fn) = [f.(fn)(:)', {fv}];
        end
    end
end


% flatten log structure
function log = flatlog(log, level)

% decompose
logc = cell(numel(log), 1);
for cc = 1:numel(log)
    logc{cc} = log(cc);
end

% for each cell, unlevel
levfield = sprintf('Level%d', level);
for cc = 1:numel(logc)

    % contains a sub-level
    if isfield(logc{cc}, levfield)

        % flatten first
        logc{cc}.(levfield) = flatlog(logc{cc}.(levfield), level + 1);

        % get fields of top struct
        topfields = fieldnames(logc{cc});
        topfields(strcmp(topfields, levfield)) = [];

        % then join
        logc{cc} = catstruct(logc{cc}, logc{cc}.(levfield));

        % then remove level
        logc{cc} = rmfield(logc{cc}, levfield);

        % and forward into lower fields
        for lcc = 2:numel(logc{cc})
            for lfc = 1:numel(topfields)
                logc{cc}(lcc).(topfields{lfc}) = logc{cc}(1).(topfields{lfc});
            end
        end
    end
end

% then join again
if numel(logc) > 1
    log = catstruct(logc{:});
else
    log = logc{1};
end
