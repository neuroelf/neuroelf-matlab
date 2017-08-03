function bdata = readbesa(filename)
% readbesa  - read in BESA export files for further processing
%
% FORMAT:       bdata = readbesa(filename)
%
% Input fields:
%
%       filename    filename of BESA export file, supported types:
%                   - avr (average file)
%                   - elp
%                   - evt (event file)
%                   - generic (control file for arbitrary dat file)
%                   - mul (multiplexed data)
%                   - pos (position file)
%                   - sfp (surface fiducial points)
%
% Output fields:
%
%       bdata       BESA file as 1x1 struct with content with fields
%        .NrOfChannels      (avr, generic, mul, pos)
%        .NrOfSamples       (avr, generic, mul)
%        .SamplingFrequency (generic, mul, in Hz)
%        .ChannelLabels     (avr, mul)
%        .AverageData       (avr, CxT)
%        .ChannelData       (mul, CxT)
%        .GenericData       (generic, CxT, transio!)
%        .DataType          (generic, used for different datatypes)
%        .PositionData      (pos, Cx9)
%        .Channel           (elp)
%        .NrOfEvents        (evt)
%        .Event             (evt)
%        .NrOfFiducials     (sfp)
%        .Fiducial          (sfp)
%        .FirstLine         copy of first line of file
%
% Note: the coordinate system will be automatically transformed into
%       a mm-unit space!
%
% Note: text files > 16MB are not supported at this time!

% Version:  v0.9b
% Build:    11050712
% Date:     Apr-08 2011, 9:16 PM EST
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
if nargin < 1 || ...
   ~ischar(filename) || ...
    isempty(filename) || ...
    exist(filename(:)', 'file') ~= 2
    error( ...
        'neuroelf:BadArgument',...
        'Bad or missing argument. Try ''help %s''.',...
        mfilename ...
    );
end

% read in file
filename = filename(:)';
try
    fileinfo = dir(filename);
    if numel(fileinfo) ~= 1 || ...
        fileinfo.isdir
        error('NOT_A_VALID_FILE');
    end
    if fileinfo.bytes > 16777216
        error('FILE_TOO_LARGE');
    end
    filecont = splittocell(asciiread(filename), char([10, 13]), 1, 1);
catch ne_eo;
    error( ...
        'neuroelf:FileReadError', ...
        'Error reading file: ''%s'' (%s).', ...
        filename, ne_eo.message ...
    );
end

% get first line to determine format
firstline = filecont{1};
if ~isempty(regexpi(firstline, 'nchan=\s*\d+\s+segmentname='));
    itype = 'avr';
elseif ~isempty(regexpi(firstline, '^tmu\s+code\s+trino'))
    itype = 'evt';
elseif ~isempty(regexpi(firstline, 'besa\s+generic\s+data'))
    itype = 'generic';
elseif ~isempty(regexpi(firstline, ...
    '^timepoints=\s*\d+\s+channels=\s*\d+\s+beginsweep'))
    itype = 'mul';
elseif ~isempty(regexpi(firstline, ...
    '^[a-z_][a-z_0-9]*\s+\-?0\.[0-9]+\s+\-?0\.[0-9]+\s+\-?0\.[0-9]+\s*$'))
    itype = 'sfp';
elseif ~isempty(regexpi(firstline, ...
    '^\s*(\-?[01]\.[0-9]+\s+){2,8}(\-?[01]\.[0-9]+)\s*$')) || ...
   ~isempty(regexpi(firstline, ...
    '^\s*(\-?[0-9]+\.[0-9]+\s+){2,8}(\-?[0-9]\.[0-9]+)\s*$'))
    itype = 'pos';
elseif ~isempty(regexpi(firstline, '^(eeg|icr|meg|pol)\s+'))
    itype = 'elp';
else
    error( ...
        'neuroelf:BadArgument', ...
        'The specified file couldn''t be identified.' ...
    );
end

% initialize output
bdata = struct;

% reading depends on type
switch (itype)

    % AVR files
    case {'avr'}

        % get settings
        [rxm{1:3}] = regexpi(firstline, 'npts=\s*(\d+)\s+');
        if isempty(rxm{3})
            error( ...
                'neuroelf:BadFileContent', ...
                'Couldn''t get number of points from file.' ...
            );
        end
        try
            numsamples = str2double(firstline(rxm{3}{1}(1):rxm{3}{1}(2)));
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            error( ...
                'neuroelf:BadFileContent', ...
                'Couldn''t convert number of points from file.' ...
            );
        end
        [rxm{1:3}] = regexpi(firstline, 'nchan=\s*(\d+)\s+');
        try
            numchannels = str2double(firstline(rxm{3}{1}(1):rxm{3}{1}(2)));
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            error( ...
                'neuroelf:BadFileContent', ...
                'Couldn''t convert number of channels from file.' ...
            );
        end
        [rxm{1:3}] = regexpi(firstline, '\s+di=\s*([0-9]+\.?[0-9]*)\s+');
        try
            sfreq = 1000 / str2double(firstline(rxm{3}{1}(1):rxm{3}{1}(2)));
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            error( ...
                'neuroelf:BadFileContent', ...
                'Couldn''t convert sampling frequency from file.' ...
            );
        end
        bdata.NrOfChannels = numchannels;
        bdata.NrOfSamples = numsamples;
        bdata.SamplingFrequency = sfreq;
        bdata.FirstLine = firstline;

        % check file size
        if numel(filecont) < (numchannels + 2)
            error( ...
                'neuroelf:BadFileContent', ...
                'Too few lines in file to read data.' ...
            );
        end

        % read channel names
        clabels = splittocell(filecont{2}, ' ');
        if isempty(clabels{1})
            clabels(1) = [];
        end
        if numel(clabels) ~= numchannels
            error( ...
                'neuroelf:BadFileContent', ...
                'Incorrect number of channel labels in file.' ...
            );
        end
        bdata.ChannelLabels = clabels(:);

        % read channel data
        cdata = zeros(numchannels, numsamples);
        for lc = 1:numchannels
            try
                eval(['cdata(lc, :) = [' filecont{lc + 2} '];']);
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
                error( ...
                    'neuroelf:BadFileContent', ...
                    'Error converting samples for channel %d.', ...
                    lc ...
                );
            end
        end
        bdata.AverageData = cdata;


    % ELP files
    case {'elp'}

        % read every line
        cc = 0;
        bdata.NrOfChannels = cc;
        for lc = 1:numel(filecont)
            channelconf = splittocell(filecont{lc}, char([9, 32]), 1, 1);
            chconf = struct;
            numfields = numel(channelconf);
            if numfields > 0
                chconf.ChannelType = channelconf{1};
            else
                continue;
            end
            if numfields > 1
                chconf.ChannelLabel = channelconf{2};
            else
                chconf.ChannelLabel = sprintf('Channel_%d', lc);
            end
            if numfields > 3 && ...
               ~isempty(regexpi(channelconf{3}, '[\+\-]?[0-9\.e\+\_]+')) && ...
               ~isempty(regexpi(channelconf{4}, '[\+\-]?[0-9\.e\+\_]+'))
                chconf.ChannelValues = ...
                    [str2double(channelconf{3}), str2double(channelconf{4})];
            else
                chconf.ChannelValues = [];
            end
            cc = cc + 1;
            bdata.Channel(cc) = chconf;
        end
        bdata.NrOfChannels = cc;
        clabels = cell(cc, 1);
        for lc = 1:cc
            clabels{lc} = bdata.Channel(lc).ChannelLabel;
        end
        bdata.ChannelLabels = clabels;


    % EVT files
    case {'evt'}

        % read events for lines 2 to end
        ec = 0;
        bdata.NrOfEvents = ec;
        for lc = 2:numel(filecont)
            [rxm{1:3}] = regexpi(filecont{lc}, ...
                '^\s*(\d+)\s+(\d+)\s+(\d+)\s+([a-z0-9_].*)$');
            if ~isempty(rxm{3}) && ...
               all(size(rxm{3}{1}) == [4, 2])
                event = struct;
                event.TimeMu = ...
                    str2double(filecont{lc}(rxm{3}{1}(1, 1):rxm{3}{1}(1, 2)));
                event.EventCode = ...
                    str2double(filecont{lc}(rxm{3}{1}(2, 1):rxm{3}{1}(2, 2)));
                event.TrialNumber = ...
                    str2double(filecont{lc}(rxm{3}{1}(3, 1):rxm{3}{1}(3, 2)));
                event.Comment = ...
                    deblank(filecont{lc}(rxm{3}{1}(4, 1):rxm{3}{1}(4, 2)));
                ec = ec + 1;
                bdata.Event(ec) = event;
            end
        end
        bdata.NrOfEvents = ec;


    % Generic files
    case {'generic'}

        % we might need a path, so get that first
        fpath = fileparts(filename);
        if isempty(fpath)
            fpath = '.';
        end

        % read fields of lines 2...N
        bdata.NrOfChannels = 0;
        bdata.NrOfSamples = 0;
        bdata.SamplingFrequency = 1;
        bdata.Datatype = 'UNKNOWN';
        bdata.SourceFile = '';
        for lc = 2:numel(filecont)
            linec = filecont{lc};
            eqs = find(linec == '=');
            if isempty(eqs) || ...
                eqs(1) == 1
                continue;
            end
            fieldn = linec(1:(eqs(1) - 1));
            switch (lower(fieldn))
                case {'file'}
                    bdata.SourceFile = [fpath '/' linec((eqs(1) + 1):end)];
                case {'format'}
                    bdata.Datatype = linec((eqs(1) + 1):end);
                case {'nchannels'}
                    bdata.NrOfChannels = str2double(linec((eqs(1) + 1):end));
                case {'nsamples'}
                    bdata.NrOfSamples = str2double(linec((eqs(1) + 1):end));
                case {'srate'}
                    bdata.SamplingFrequency = ...
                        str2double(linec((eqs(1) + 1):end));
                otherwise
                    if strcmp(makelabel(fieldn), fieldn)
                        bdata.(fieldn) = linec((eqs(1) + 1):end);
                    end
            end
        end

        % check some fields
        switch (lower(bdata.Datatype))
            case {'float', 'float32', 'single'}
                bdata.Datatype = 'single';
            otherwise
                error( ...
                    'neuroelf:BadFileContent', ...
                    'Error reading or invalid datatype.' ...
                );
        end
        if bdata.NrOfChannels == 0 || ...
            bdata.NrOfSamples == 0
            error( ...
                'neuroelf:BadFileContent', ...
                'Error getting number of channels/samples.' ...
            );
        end

        % create transio object for data access
        try
            if ~isempty(bdata.SourceFile)
                bdata.GenericData = ...
                    transio(bdata.SourceFile, 'ieee-le', 'single', 0, ...
                        [bdata.NrOfChannels, bdata.NrOfSamples]);
            else
                bdata.GenericData = 'NO_FILE_GIVEN';
            end
        catch ne_eo;
            bdata.GenericData = sprintf('transio error: %s', ne_eo.message);
            error( ...
                'neuroelf:BadLinkedFile', ...
                'Improper linked file: ''%s''; error: %s.', ...
                bdata.SourceFile, ne_eo.message ...
            );
        end


    % MUL files
    case {'mul'}

        % read settings
        bdata.NrOfChannels = 0;
        bdata.NrOfSamples = 0;
        bdata.SamplingFrequency = 1;
        bdata.Time = '00:00:00';
        bdata.FirstLine = firstline;
        bdata.ChannelLabels = {};

        % parse settings
        settings = splittocell(firstline, char([9, 32]), 1, 1);
        if isempty(settings{1})
            settings(1) = [];
        end
        for sc = 1:numel(settings)
            setting = splittocell(settings{sc}, '=');
            if numel(setting) == 2 && ...
               ~isempty(setting{1})
                switch (lower(setting{1}))
                    case {'channels'}
                        bdata.NrOfChannels = str2double(setting{2});
                    case {'samplinginterval[ms]'}
                        bdata.SamplingFrequency = 1000 / str2double(setting{2});
                    case {'timepoints'}
                        bdata.NrOfSamples = str2double(setting{2});
                    otherwise
                        bdata.(makelabel(setting{1})) = setting{2};
                end
            end
        end

        % check
        if bdata.NrOfChannels == 0 || ...
            bdata.NrOfSamples == 0
            error( ...
                'neuroelf:BadFileContent', ...
                'Error parsing number of channels/samples.' ...
            );
        end

        % check file size
        if numel(filecont) < (bdata.NrOfSamples + 2)
            error( ...
                'neuroelf:BadFileContent', ...
                'Too few lines to read data samples.' ...
            );
        end
        cdata = zeros(bdata.NrOfSamples, bdata.NrOfChannels);

        % read channel labels
        clabels = splittocell(filecont{2}, char([9, 32]), 1, 1);
        if isempty(clabels{1})
            clabels(1) = [];
        end
        if numel(clabels) ~= bdata.NrOfChannels
            error( ...
                'neuroelf:BadFileContent', ...
                'Incorrect number of channel labels in file.' ...
            );
        end
        bdata.ChannelLabels = clabels(:);

        % read channel data
        for lc = 1:bdata.NrOfSamples
            try
                eval(['cdata(lc, :) = [' filecont{lc + 2} '];']);
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
                error( ...
                    'neuroelf:BadFileContent', ...
                    'Error converting data for sample %d.', ...
                    lc ...
                );
            end
        end
        bdata.ChannelData = cdata';


    % POS files
    case {'pos'}

        % only very few fields...
        bdata.NrOfChannels = numel(filecont);
        posdata = zeros(numel(filecont), 9);

        % reading pos data
        for lc = 1:numel(filecont)
            try
                eval(['posline = [' filecont{lc} '];']);
                posdata(lc, 1:min(9, numel(posline))) = posline;
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
                error( ...
                    'neuroelf:BadFileContent', ...
                    'Invalid position data in line %d.', ...
                    lc ...
                );
            end
        end
        if all(lsqueeze(abs(posdata(:, 1:6))) < 2)
            posdata(:, 1:6) = 1000 * posdata(:, 1:6);
        end
        bdata.PositionData = posdata;


    % SFP files
    case {'sfp'}

        % only very few fields...
        bdata.NrOfFiducials = numel(filecont);

        % parse fiducials
        fc = 0;
        for lc = 1:numel(filecont)
            linec = filecont{lc};
            spp = find(linec == char(9) | linec == char(32));
            if ~isempty(spp) && ...
                spp(1) > 1 && ...
                strcmp(makelabel(linec(1:(spp(1) - 1))), linec(1:(spp(1) - 1)))
                fid.Label = linec(1:(spp(1) - 1));
                try
                    eval(['fid.Coordinate = 1000 * [' linec((spp(1) + 1):end) '];']);
                    fc = fc + 1;
                    bdata.Fiducial(fc) = fid;
                catch ne_eo;
                    neuroelf_lasterr(ne_eo);
                end
            end
        end
        bdata.NrOfFiducials = fc;


    otherwise
        error( ...
            'neuroelf:ImplementationError', ...
            'Error in M-file %s.', ...
            mfilename ...
        );
end
