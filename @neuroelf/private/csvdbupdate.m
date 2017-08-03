function csvdbupdate(dbfile, inopts)
% csvdbupdate  - update the data field from CSV files
%
% FORMAT:       csvdbupdate(dbfile, opts)
%
% Input fields:
%
%       dbfile      filename of MAT file containing the DB
%       opts        optional settings (overriding opts in DB file)
%        .convert   convert flag for acsvread
%        .csvinpat  filename pattern for input CSV files
%        .outfile   write text based output to file
%        .spssout   boolean flag whether or not to write SAV output
%        .validIDs  cell array with IDs (to add)
%        .writeout  boolean flag whether or not to write output
%        .xlsout    boolean flag whether or not to write XLS output
%
% No output fields.

% Version:  v0.9b
% Build:    11050712
% Date:     Apr-09 2011, 1:55 PM EST
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

% argument check
if nargin < 1 || ...
   ~ischar(dbfile) || ...
    isempty(dbfile) || ...
    exist(dbfile(:)', 'file') ~= 2
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing dbfile argument.' ...
    );
end

% try loading DB
opwd = pwd;
dbfile = dbfile(:)';
try
    db = load(dbfile);
    if any(dbfile == filesep)
        try
            cd(fileparts(dbfile));
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            warning( ...
                'neuroelf:CDError', ...
                'Error changing directory to ''%s''.', ...
                fileparts(dbfile) ...
            );
        end
    end
    data = db.data;
    opts = db.opts;
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    error( ...
        'neuroelf:DBLoadError', ...
        'Error loadind CSV based DB from file.' ...
    );
end

% update options
if nargin < 2 || ...
   ~isstruct(inopts) || ...
    numel(inopts) ~= 1
    inopts = struct;
end
ofc = fieldnames(opts);
for fc = 1:numel(ofc)
    if isfield(inopts, ofc{fc}) && ...
        strcmp(class(inopts.(ofc{fc})), class(opts.(ofc{fc})))
        if iscell(opts.(ofc{fc}))
            if isnumeric(inopts.(ofc{fc}){1})
                opts.(ofc{fc}) = {union([opts.(ofc{fc}){:}], [inopts.(ofc{fc}){:}])};
            else
                opts.(ofc{fc}) = union(opts.(ofc{fc})(:), inopts.(ofc{fc})(:));
            end
        else
            opts.(ofc{fc}) = inopts.(ofc{fc})(:)';
        end
    end
end

% get ID column in data
didc = find(strcmp(data(1, :)', 'ID'));
if numel(didc) ~= 1
    error( ...
        'neuroelf:DBError', ...
        'Corrupt data table, no ID column found.' ...
    );
end

% work on list
nd = {''};
nd = nd(1, ones(1, size(data, 2)));
if any(opts.csvinpat == filesep)
    cfld = fileparts(opts.csvinpat);
    cpat = strrep(opts.csvinpat, [cfld filesep], '');
else
    cfld = pwd;
    cpat = opts.csvinpat;
end
csvs = findfiles(cfld, cpat, 'depth=1');
for fc = 1:numel(csvs)

    % try loading data
    try
        if opts.convert
            csv = acsvread(csvs{fc}, '","', struct('convert', 1));
        else
            csv = acsvread(csvs{fc}, '","');
        end
        csva = splittocell(csv{1}{1}, ',');
        csva{end} = strrep(csva{end}, '"', '');
        csv{1} = [csva, csv{1}(2:end+1-numel(csva))];
        for lc = 1:numel(csv)
            csv{lc}{1} = regexprep(csv{lc}{1}, '^\"', '');
            csv{lc}{end} = regexprep(csv{lc}{end}, '\",?$', '');
        end

        tt = regexprep(csv{1}(:), '^.*<(.*)>$', '$1');
        idc = find(strcmp(tt, 'ID'));
        if numel(idc) ~= 1
            error('No ID column.');
        end
    catch ne_eo;
        warning( ...
            'neuroelf:DBFileError', ...
            'Error reading file %s: %s.', ...
            csvs{fc}, ne_eo.message ...
        );
        continue;
    end

    % create overlap indices
    [itag, ita, itb] = intersect(data(1, :)', tt);
    if numel(itag) ~= numel(tt)
        warning( ...
            'neuroelf:DBFileError', ...
            'Superfluous/double tags found in ''%s'', will be dropped.', ...
            csvs{fc} ...
        );
    end

    % go through data
    for lc = 2:numel(csv)

        % convert numerical data
        csvl = regexprep(csv{lc}(:), ...
            '^\s*([+\-]?\d+(\.\d+)?([eE][+\-]?\d+)?)\s*$', '==$1==');
        for cc = 1:numel(csvl)
            if numel(csvl{cc}) > 4 && ...
                all(csvl{cc}([1:2, end-1:end]) == '=')
                csvl{cc} = str2double(csvl{cc}(3:end-2));
            end
        end
        csv{lc}(:) = csvl(:);

        % get ID
        id = csv{lc}{idc};
        if isempty(id)
            continue;
        end
        idf = false;
        for ic = 1:numel(opts.validIDs{1})
            if isequal(opts.validIDs{1}(ic), id)
                idf = true;
                break;
            end
        end
        if ~idf
            disp(sprintf('ID %d not in list of valid IDs.', id));
            continue;
        end

        % try to find corresponding row in dataset
        idf = false;
        for rc = 2:size(data, 1)
            if isequal(id, data{rc, didc})
                idf = true;
                break;
            end
        end

        % if found
        if idf

            % update row
            data(rc, ita) = csv{lc}(itb);

        % otherwise
        else

            % create new row
            data(end + 1, :) = nd;
            data(end, ita) = csv{lc}(itb);
        end
    end
end

% save updated data
try
    cd(opwd);
catch ne_eo;
    neuroelf_lasterr(ne_eo);
end
try
    save(dbfile, 'data', 'opts', '-v6');
catch ne_eo;
    warning( ...
        'neuroelf:DBSaveError', ...
        'Error saving DB back to file %s (%s).', ...
        dbfile, ne_eo.message ...
    );
    return;
end

% if output requested, try output
if ~opts.writeout
    return;
end
try
    cd(fileparts(dbfile));
catch ne_eo;
    neuroelf_lasterr(ne_eo);
end
try
    fid = fopen(opts.outfile, 'w');
    if fid < 1
        error('File not writable.');
    end

    % write header
    disp('Saving CSV file...');
    tagstring = sprintf('"%s"\t', data{1, :});
    fprintf(fid, '%s\n', tagstring(1:end-1));

    % write contents of DB
    for rc = 2:size(data, 1)
        dataline = {''};
        dataline = dataline(1, ones(1, size(data, 2)));
        for ec = 1:size(data, 2)
            if ischar(data{rc, ec})
                if regexp(data{rc, ec}, '^[+\-]?\d+(\.\d+)?([eE][+\-]\d+)?$')
                    dataline{ec} = data{rc, ec};
                elseif ~isempty(data{rc, ec})
                    dataline{ec} = ['"' data{rc, ec} '"'];
                end
            else
                dataline{ec} = sprintf('%d', data{rc, ec});
            end
        end
        datastring = sprintf('%s\t', dataline{:});
        fprintf(fid, '%s\n', datastring(1:end-1));
    end

    % close file
    fclose(fid);
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    error( ...
        'neuroelf:FileOpenError', ...
        'Output file not writable.' ...
    );
end

% also write SPSS/SAV output
if opts.spssout
    try
        s = [];
        s = xff('new:spss');
        s.ImportData(data(2:end, :), struct('vars', {data(1, :)}));
        disp('Saving SPSS file...');
        s.SaveAs([opts.outfile(1:end-3) 'sav']);
        s.ClearObject;
    catch ne_eo;
        oeo = ne_eo;
        if isxff(s)
            try
                s.ClearObject;
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
            end
        end
        warning( ...
            'neuroelf:DBError', ...
            'Error writing XLS output sheet (%s).', ...
            oeo.message ...
        );
    end
end

% also write XLS output
if opts.xlsout
    try
        disp('Saving XLS file...');
        xlswrite([opts.outfile(1:end-3) 'xls'], data');
    catch ne_eo;
        warning( ...
            'neuroelf:DBError', ...
            'Error writing XLS output sheet (%s).', ...
            ne_eo.message ...
        );
    end
end

% go back to old folder
try
    cd(opwd);
catch ne_eo;
    neuroelf_lasterr(ne_eo);
end
