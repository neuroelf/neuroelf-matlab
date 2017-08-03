function xo = spss_ImportData(xo, data, opts)
% SPSS::ImportData  - import data into an SPSS file object
%
% FORMAT:       [spss] = spss.ImportData(data [, opts])
%
% Input fields:
%
%       data        CxV cell array with C cases and V variables
%       opts        optional settings
%        .defpform  default PrintFormat for numeric variables (329731)
%        .ovwrite   boolean flag to overwrite existing data, default: true
%        .vars      1xV variable names, if not given {'VAR00001', ...}
%        .varlabels 1xV variable labels (auto-determined for 'AUTO' type)
%                   cells must be a Lx2 array with {value, 'Label'; ...}
%        .vartypes  1xV variable types, if not given {'AUTO', ...}
%                   selection of {'AUTO'}, 'NUMERIC', 'STRING'
%
% Output field:
%
%       spss        altered object
%
% Using: gluetostring, makelabel.

% Version:  v1.1
% Build:    16020314
% Date:     Feb-03 2016, 2:44 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2014, 2016, Jochen Weber
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

% neuroelf library
global ne_methods;

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'spss') || ...
   (~iscell(data) && ~isa(data, 'double')) || isempty(data) || ndims(data) > 2
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
if ~iscell(data)
    dcell = cell(size(data));
    for dcc = 1:numel(dcell)
        dcell{dcc} = data(dcc);
    end
    data = dcell;
end
bc = xo.C;
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'defpform') || ~isa(opts.defpform, 'double') || numel(opts.defpform) ~= 1 || ...
    isnan(opts.defpform) || isinf(opts.defpform) || ...
    opts.defpform ~= round(opts.defpform) || opts.defpform < 0 || ...
    any(sprintf('%08x', opts.defpform) < '00010101') || ...
    any(sprintf('%08x', opts.defpform) > '001f1f1f')
    opts.defpform = 329731;
end
if ~isfield(opts, 'ovwrite') || ~islogical(opts.ovwrite) || numel(opts.ovwrite) ~= 1
    opts.ovwrite = true;
end
if ~isfield(opts, 'vars') || ~iscell(opts.vars) || numel(opts.vars) ~= size(data, 2)
    opts.vars = cell(1, size(data, 2));
end
if ~isfield(opts, 'varlabels') || ~iscell(opts.varlabels) || numel(opts.varlabels) ~= size(data, 2)
    opts.varlabels = cell(1, size(data, 2));
else
    opts.varlabels = opts.varlabels(:)';
end
if ~isfield(opts, 'varnames') || ~iscell(opts.varnames) || ...
   ~isequal(size(opts.vars), size(opts.varnames))
    opts.varl = cell(1, numel(opts.vars));
else
    opts.varl = opts.varnames(:)';
end
opts.vars = opts.vars(:)';
makelabel = ne_methods.makelabel;
for vc = 1:numel(opts.vars)
    if ~ischar(opts.vars{vc}) || isempty(opts.vars{vc})
        opts.vars{vc} = sprintf('VAR%05d', vc);
    elseif ~strcmp(makelabel(opts.vars{vc}(:)'), opts.vars{vc}(:)')
        opts.vars{vc} = makelabel(opts.vars{vc}(:)');
    end
    if numel(opts.vars{vc}) > 8 || ~all(opts.vars{vc} == upper(opts.vars{vc}))
        opts.varl{vc} = opts.vars{vc};
        opts.vars{vc} = upper(opts.vars{vc}(1:min(8, numel(opts.vars{vc}))));
    end
end
uvars = unique(opts.vars(:));
if numel(uvars) ~= numel(opts.vars)
    error('neuroelf:xff:badArgument', 'Incongruent variable naming.');
end
if ~isfield(opts, 'vartypes') || ~iscell(opts.vartypes) || numel(opts.vartypes) ~= size(data, 2)
    opts.vartypes = {'AUTO'};
    opts.vartypes = opts.vartypes(1, ones(1, size(data, 2)));
end
opts.vartypes = opts.vartypes(:)';
for vc = 1:numel(opts.vartypes)
    if ~ischar(opts.vartypes{vc}) || ~any(strcmpi(opts.vartypes{vc}(:)', {'AUTO', 'NUMERIC', 'STRING'}))
        opts.vartypes{vc} = 'AUTO';
    else
        opts.vartypes{vc} = upper(opts.vartypes{vc}(:)');
    end
end

% determine auto types
sar = {''};
sar = sar(ones(1, size(data, 1)), 1);
for vc = 1:numel(opts.vartypes)
    if strcmp(opts.vartypes{vc}, 'AUTO')

        % look at all data and update numeric/string counters
        nc = 0;
        sc = 0;
        sa = sar;
        for rc = 1:size(data, 1)
            if ischar(data{rc, vc}) && ~isempty(data{rc, vc})
                sa{rc} = data{rc, vc}(:)';
                sc = sc + 1;
            elseif isa(data{rc, vc}, 'double') && numel(data{rc, vc}) == 1
                nc = nc + 1;
            end
        end

        % make decision for numeric if at least as many numbers as strings
        if nc >= sc
            opts.vartypes{vc} = 'NUMERIC';

            % make sure everything is a 1x1 double
            for rc = 1:size(data, 1)
                if ~isa(data{rc, vc}, 'double') || isempty(data{rc, vc})
                    data{rc, vc} = NaN;
                else
                    data{rc, vc} = data{rc, vc}(1);
                end
            end

        % otherwise
        else

            % make sure everything is a string first
            for rc = 1:size(data, 1)
                if ~ischar(data{rc, vc}) || isempty(data{rc, vc})
                    data{rc, vc} = '';
                else
                    data{rc, vc} = data{rc, vc}(:)';
                end
            end

            % check for regular labels
            usa = unique(sa);
            usc = numel(usa);
            if usc <= sqrt(sc)
                opts.vartypes{vc} = 'NUMERIC';

                % for only two labels other than '' use 0, 1
                if sum(~strcmp(usa, '')) == 2
                    lcc = 0;

                % multiple strings > labels from 1...N
                else
                    lcc = 1;
                end

                % find labels and replace in data
                labidx = cell(1, usc);
                for lc = 1:usc
                    labidx{lc} = find(strcmp(data(:, vc), usa{lc}));
                end
                labval = cell(0, 2);
                for lc = 1:usc
                    if isempty(usa{lc})
                        data(labidx{lc}, vc) = repmat({NaN}, numel(labidx), 1);
                    else
                        data(labidx{lc}, vc) = repmat({lcc}, numel(labidx), 1);
                        labval(end + 1, :) = {lcc, usa{lc}};
                        lcc = lcc + 1;
                    end
                end
                opts.varlabels = labval;

            % otherwise replace non-chars and empty strings with ''
            else
                opts.vartypes{vc} = 'STRING';
            end
        end
    end

    % detect required length
    if strcmp(opts.vartypes{vc}, 'STRING')
        opts.varlabels{vc} = 8 * ceil(size(char(data(:, vc)), 2) / 8);
    end
end

% overwrite old data?
if opts.ovwrite

    % set fields in FileHeader
    bc.FileHeader.NrOfCases = size(data, 1);

    % overwrite variables
    bc.Variables(:) = [];
    bc.ValueLabels(:) = [];
    for rc = 1:numel(bc.InfoRecords)
        if bc.InfoRecords(rc).SubType == 11
            bc.InfoRecords(rc).NrOfElements = 0;
            bc.InfoRecords(rc).Data = uint8(zeros(4, 0));
            bc.InfoRecords(rc).Parsed(:) = [];
        elseif bc.InfoRecords(rc).SubType == 13
            bc.InfoRecords(rc).NrOfElements = 0;
            bc.InfoRecords(rc).Data = uint8(zeros(1, 0));
            bc.InfoRecords(rc).Parsed = struct;
        elseif bc.InfoRecords(rc).SubType == 16
            bc.InfoRecords(rc).NrOfElements = 2;
            bc.InfoRecords(rc).Data = uint8(zeros(8, 2));
            bc.InfoRecords(rc).Parsed = struct;
        end
    end

    % put data into data
    bc.DataRecords = data;

    % count unique labels
    vct = numel(opts.vars);
    ulab = {};
    vlab = {};
    for vc = 1:vct
        if ~isempty(opts.varlabels{vc}) && iscell(opts.varlabels{vc})
            ulaba = any2ascii(opts.varlabels{vc});
            ulabp = find(strcmp(ulab, ulaba));
            if isempty(ulabp)
                ulab{end + 1} = ulaba;
                vlab{end + 1} = opts.varlabels{vc};
                opts.varlabels{vc} = -numel(ulab);
            else
                opts.varlabels{vc} = -ulabp(1);
            end
        end
    end

    % re-fill Variables
    bc.Variables(vct).RecordType = 2;
    rst11d = repmat(uint8([3, 8, 1; zeros(3, 3)]), 1, vct);
    rst11p = struct('VarMeasure', 1, 'VarWidth', 8, 'VarAlignment', 1);
    rst11p = rst11p(1, ones(1, vct));
    for vc = 1:vct
        bc.Variables(vc).RecordType = 2;
        if strcmp(opts.vartypes{vc}, 'NUMERIC')
            bc.Variables(vc).VarType = 0;
        else
            bc.Variables(vc).VarType = opts.varlabels{vc};
            rst11d(1, vc * 3 - 1) = opts.varlabels{vc};
            rst11p(vc).VarWidth = opts.varlabels{vc};
        end
        bc.Variables(vc).HasLabel = 0;
        bc.Variables(vc).NrOfMissingValues = 0;
        bc.Variables(vc).PrintFormat = opts.defpform;
        bc.Variables(vc).WriteFormat = opts.defpform;
        bc.Variables(vc).Name = opts.vars{vc};
    end
    for rc = 1:numel(bc.InfoRecords)
        if bc.InfoRecords(rc).SubType == 11
            bc.InfoRecords(rc).NrOfElements = size(rst11d, 2);
            bc.InfoRecords(rc).Data = uint8(rst11d);
            bc.InfoRecords(rc).Parsed = rst11p;
        elseif bc.InfoRecords(rc).SubType == 16
            bc.InfoRecords(rc).Data(1, 1) = 1;
            nocas = size(data, 1);
            bc.InfoRecords(rc).Data(1, 2) = mod(nocas, 256);
            bc.InfoRecords(rc).Data(2, 2) = mod(floor(nocas / 256), 256);
            bc.InfoRecords(rc).Data(3, 2) = mod(floor(nocas / 65536), 256);
        end
    end

    % re-fill ValueLabels
    if ~isempty(ulab)
        bc.ValueLabels(numel(ulab)).RecordType = 3;
    end
    for vc = 1:numel(ulab)
        bc.ValueLabels(vc).RecordType = 3;
        bc.ValueLabels(vc).Labels = struct('Value', 0, 'LabelLen', 0, 'Label', '');
        for lc = 1:size(vlab{vc}, 1)
            bc.ValueLabels(vc).Labels(lc).Value = vlab{vc}{lc, 1};
            bc.ValueLabels(vc).Labels(lc).Label = vlab{vc}{lc, 2};
            bc.ValueLabels(vc).Labels(lc).LabelLen = numel(vlab{vc}{lc, 2});
        end
        vmat = false(1, vct);
        for vcc = 1:vct
            if opts.varlabels{vcc} == -vc
                vmat(vcc) = true;
            end
        end
        bc.ValueLabels(vc).Variables = find(vmat);
    end

    % create longer labels array
    llab = struct;
    for vc = 1:vct
        if ~isempty(opts.varl{vc})
            llab.(opts.vars{vc}) = opts.varl{vc};
        end
    end
    llabf = fieldnames(llab);
    for vc = 1:numel(llabf)
        llabf{vc} = [llabf{vc} '=' llab.(llabf{vc})];
    end
    if ~isempty(llabf)
        llabf = ne_methods.gluetostring(llabf, char(9));
        llabp = mod(numel(llabf), 8);
        if llabp > 0
            llabf(end+1:end+8-llabp) = ' ';
        end
        for rc = 1:numel(bc.InfoRecords)
            if bc.InfoRecords(rc).SubType == 13
                bc.InfoRecords(rc).NrOfElements = numel(llabf);
                bc.InfoRecords(rc).Data = uint8(double(llabf));
                bc.InfoRecords(rc).Parsed = llab;
                break;
            end
        end
    end

% try to add to data
else

    % for now only give a warning and return
    warning('neuroelf:xff:notYetImplemented', ...
        'Appending of imported data not yet implemented.');
    return;
end

% set in memory
xo.C = bc;
