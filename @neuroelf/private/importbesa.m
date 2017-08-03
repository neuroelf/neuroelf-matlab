function [cdd, ccc, ctc, act, prt] = importbesa(pattern)
% importbesa  - import EEG/MEG data from a set of BESA files
%
% FORMAT:       [cdd, ccc, ctc, act, prt] = importbesa(pattern)
%
% Input fields:
%
%       pattern     BESA exported files pattern (e.g. 'run1.*')
%
% Output fields:
%
%       cdd         channel data defition object
%       ccc         channel coordinate configuration object
%       ctc         channel time course object
%       act         average channel time course object
%       prt         protocol (events) object
%
% See also readbesa

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:20 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2016, Jochen Weber
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
   ~ischar(pattern) || ...
    isempty(pattern)
    error( ...
        'neuroelf:BadArgument',...
        'Bad or missing argument. Try ''help %s''.',...
        mfilename ...
    );
end

% find matching files
[ppath{1:3}] = fileparts(pattern(:)');
if isempty(ppath{1})
    ppath{1} = '.';
end
besafiles = findfiles(ppath{1}, [ppath{2} '.*']);
if numel(besafiles) < 1
    error( ...
        'neuroelf:FileNotFound', ...
        'No matching files found.' ...
    );
end

% create objects
cdd = xff('new:cdd');
ccc = xff('new:ccc');
ctc = xff('new:ctc');
act = xff('new:ctc');
prt = xff('new:prt');
ols = {cdd, ccc, ctc, act, prt};

% tctype predefinitions
tctraw = false;
tctavg = false;
hasprt = false;

% try to parse files
readfields = struct;
for fc = 1:numel(besafiles)
    try
        besafile = readbesa(besafiles{fc});
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        continue;
    end
    fields = fieldnames(besafile);
    for cc = 1:numel(fields)
        if ~isfield(readfields, fields{cc})
            readfields.(fields{cc}) = besafile.(fields{cc});
        else
            switch (lower(fields{cc}))
                case {'channellabels'}
                    if numel(readfields.(fields{cc})) ~= numel(readfields.(fields{cc})) || ...
                       ~all(strcmpi(readfields.(fields{cc}), readfields.(fields{cc})))
                        clearxffobjects(ols);
                        error( ...
                            'neuroelf:BadFileContent', ...
                            'Label(s) of channel mismatch.' ...
                        );
                    end
                case {'nrofchannels'}
                    if readfields.(fields{cc}) ~= besafile.(fields{cc})
                        if readfields.(fields{cc}) < besafile.(fields{cc}) && ...
                            isfield(readfields, 'PositionData')
                            readfields.(fields{cc}) = besafile.(fields{cc});
                        elseif readfields.(fields{cc}) > besafile.(fields{cc}) && ...
                            isfield(besafile, 'PositionData')
                        else
                            clearxffobjects(ols);
                            error( ...
                                'neuroelf:BadFileContent', ...
                                'NrOfChannels mismatch between files.' ...
                            );
                        end
                    end
                case {'nrofsamples', 'samplingfrequency'}
                    if readfields.(fields{cc}) ~= besafile.(fields{cc})
                        clearxffobjects(ols);
                        error( ...
                            'neuroelf:BadFileContent', ...
                            '%s mismatch between files.', ...
                            fields{cc} ...
                        );
                    end
                otherwise
                    warning( ...
                        'neuroelf:BadFileContent', ...
                        'Ambiguous field (%s) found in file ''%s''.', ...
                        fields{cc}, ...
                        besafiles{fc} ...
                    );
            end
        end
    end
end

% check whether necessary fields are in place
if ~isfield(readfields, 'NrOfChannels') || ...
   ~isfield(readfields, 'NrOfSamples') || ...
   ~isfield(readfields, 'SamplingFrequency') || ...
   ~isfield(readfields, 'ChannelLabels') || ...
   (~isfield(readfields, 'AverageData') && ...
    ~isfield(readfields, 'ChannelData') && ...
    ~isfield(readfields, 'GenericData'))
	clearxffobjects(ols);
    error( ...
        'neuroelf:MissingData', ...
        'Some crucial information is missing, aborting import.' ...
    );
end

% set crucial fields in CDD file
cdd.NrOfChannels = readfields.NrOfChannels;
cdd.NrOfSamples = readfields.NrOfSamples;
cdd.SamplingFrequency = readfields.SamplingFrequency;
cdd.TimeOffsetMu = 0;
% cdd.TimeOffsetMu = ... (to be implemented later)
cdd.ChannelLabels = cell(numel(readfields.ChannelLabels), 2);
cdd.ChannelLabels(:, 2) = readfields.ChannelLabels(:);
cdd.NrOfPrecalculatedFiles = 0;
cdd.PrecalculatedFiles = cell(0, 1);
cdd.NrOfMarkers = 0;

% set crucial fields in CCC file
ccc.NrOfDataChannels = readfields.NrOfChannels;
ccc.DataChannelType = zeros(1, readfields.NrOfChannels);
ccc.DataChannelLabel = readfields.ChannelLabels(:);
% ccc.ReferenceType = ... (to be implemented later)
ccc.ReferenceType = 1;
ccc.CoordinateOrigin = [0 0 0];
ccc.DataCoordinates = zeros(readfields.NrOfChannels, 9);
ccc.DataCoordinates(:, 7) = 1;

% unset default position data (if available)
if isfield(readfields, 'PositionData')
    ccc.DataCoordinates(1:size(readfields.PositionData, 1), :) = ...
        readfields.PositionData;
end

ccc.UseSphericalCoords = 0;
% ccc.NrOfDataChannelParams = ... (to be implemented later)
ccc.NrOfDataChannelParams = 0;
ccc.DataChannelParamType = zeros(1, 0);
ccc.DataChannelParameter = zeros(0, readfields.NrOfChannels);
ccc.NrOfTransformations = 0;
ccc.Transformations(1:end) = [];
ccc.NrOfVirtualConfigs = 0;
ccc.VirtualConfig(1:end) = [];

% if there is channel type data, put it in!
if isfield(readfields, 'Channel') && ...
    isstruct(readfields.Channel) && ...
    isfield(readfields.Channel, 'ChannelType') && ...
    isfield(readfields.Channel, 'ChannelLabel') && ...
    numel(readfields.Channel) == readfields.NrOfChannels
    for cc = 1:readfields.NrOfChannels
        ch = readfields.Channel(cc);
        ctype = ch.ChannelType;
        clabl = ch.ChannelLabel;
        labmt = find(strcmpi(clabl, cdd.ChannelLabels(:, 2)));
        if numel(labmt) ~= 1
            warning( ...
                'neuroelf:BadFileContent', ...
                'Badly referenced channel %s, assuming 1-on-1 mapping.', ...
                labmt ...
            );
            labmt = cc;
        end
        switch (lower(ctype))
            % do nothing for directly supported
            case {'eeg'}
                ntype = 1;
            case {'icr'}
                ntype = 4;
            case {'meg'}
                if all(ccc.DataCoordinates(labmt, 7:9) == 0)
                    ctype = 'MAG';
                    ntype = 2;
                else
                    ctype = 'GRD';
                    ntype = 3;
                end
            case {'pol'}
                ntype = 5;
            otherwise
                clearxffobjects(ols);
                error( ...
                    'neuroelf:BadFileContent', ...
                    'Unsupported channel type: %s.', ...
                    ctype ...
                );
        end
        cdd.ChannelLabels{labmt, 1} = upper(ctype);
        ccc.DataChannelType(labmt) = ntype;
    end

    % check all types
    badtypes = find(ccc.DataChannelType == 0);
    if ~isempty(badtypes)
        clearxffobjects(ols);
        error( ...
            'neuroelf:BadFileContent', ...
            'Type not specified for channels %s.', ...
            gluetostring(cdd.ChannelLabels(badtypes, 2)', ', ') ...
        );
    end

else
    % assume EEG without additional pos
    ccc.DataChannelType(:) = 1;
end

% set crucial fields in CTC (if applicable)
ctc.CompressionType = 0;
ctc.NrOfCompressionParams = 0;
ctc.DataType = 1;
ctc.SampleOrdering = 2;
ctc.NrOfChannels = readfields.NrOfChannels;
ctc.NrOfSamples = readfields.NrOfSamples;
ctc.SamplingFrequency = readfields.SamplingFrequency;
if isfield(readfields, 'GenericData') && ...
    (readfields.NrOfSamples / readfields.SamplingFrequency) > 5
    tctraw = true;
    % ctc.DataType = ... (to be implemented later)
    ctc.CTCData = readfields.GenericData;
elseif isfield(readfields, 'ChannelData') && ...
    (readfields.NrOfSamples / readfields.SamplingFrequency) > 5
    tctraw = true;
    ctc.CTCData = single(readfields.ChannelData);
elseif isfield(readfields, 'AverageData') && ...
    (readfields.NrOfSamples / readfields.SamplingFrequency) > 5
    tctraw = true;
    ctc.CTCData = single(readfields.AverageData);
end

% set crucial fields in ACT (if applicable)
act.CompressionType = 0;
act.NrOfCompressionParams = 0;
act.DataType = 1;
act.SampleOrdering = 4;
act.NrOfChannels = readfields.NrOfChannels;
act.NrOfSamples = readfields.NrOfSamples;
act.SamplingFrequency = readfields.SamplingFrequency;
act.NrOfConditions = 1;
act.Condition(1).NrOfOnsets = 0;
act.Condition.NrOfPreTriggerSamples = 0;
act.Condition.NrOfPostTriggerSamples = readfields.NrOfSamples;
if isfield(readfields, 'AverageData') && ...
    (readfields.NrOfSamples / readfields.SamplingFrequency) <= 5
    tctavg = true;
    % act.NrOfPreTriggerSamples = ... (to be implemented later)
    act.Condition.CTCData = single(readfields.AverageData);
elseif isfield(readfields, 'ChannelData') && ...
    (readfields.NrOfSamples / readfields.SamplingFrequency) <= 5
    tctavg = true;
    act.Condition.CTCData = single(readfields.ChannelData);
elseif isfield(readfields, 'GenericData') && ...
    (readfields.NrOfSamples / readfields.SamplingFrequency) <= 5
    tctavg = true;
    act.Condition.CTCData = single(readfields.GenericData(:, :));
end

% any good data for a PRT??
if isfield(readfields, 'NrOfEvents') && ...
    isfield(readfields, 'Event')
    hasprt = true;

    % make initial settings
    prt.ResolutionOfTime = 'u';
    prt.Experiment = ppath{2};

    % get conditions
    cnd = struct;
    mark = struct;
    mc = 0;
    evt = readfields.Event;
    for ec = 1:numel(evt)
        evtime = evt(ec).TimeMu;
        numcode = evt(ec).EventCode;
        numtrial = evt(ec).TrialNumber;
        comment = evt(ec).Comment;
        if numcode == 1
            ecode = sprintf('Code_%d__Trial_%d', numcode, numtrial);
            if isfield(cnd, ecode)
                cnd.(ecode)(end + 1) = evtime;
            else
                cnd.(ecode) = evtime;
            end
        else
            ecode = sprintf('Marker_%d_%d_%d', evtime, numcode, numtrial);
            mc = mc + 1;
            mark.MarkerName = ecode;
            mark.TemporalPosition = evtime;
            mark.Duration = 0;
            mark.RGBColor = [0 255 128];
            mark.Comment = comment;
            if mc == 1
                mrks = mark;
            else
                mrks(mc) = mark;
            end
        end
    end

    % put conditions in PRT
    cfields = fieldnames(cnd);
    for cc = 1:numel(cfields)
        prt.AddCond(cfields{cc}, repmat(cnd.(cfields{cc})(:), [1, 2]));
    end

    % put markers in CDD
    if mc > 0
        cdd.NrOfMarkers = mc;
        cdd.Markers = mrks;
    end
end

% any good TC data
if ~tctraw && ...
   ~tctavg
    clearxffobjects(ols);
    error( ...
        'neuroelf:BadArguments', ...
        'No time course data among given files.' ...
    );
end

% set links
if tctraw
    cdd.CTCFile = [ppath{2} '.ctc'];
else
    cdd.CTCFile = '<none>';
end
if tctavg
    cdd.ACTFile = [ppath{2} '.act'];
else
    cdd.ACTFile = '<none>';
end
cdd.CCCFile = [ppath{2} '.ccc'];
if hasprt
    cdd.ProtocolFile = [ppath{2} '.prt'];
end

% save imported files
try
    cdd.SaveAs([ppath{1} '/' ppath{2} '.cdd']);
catch ne_eo;
    warning( ...
        'neuroelf:ErrorSavingXff', ...
        'Error saving xff of type CDD (%s).', ...
        ne_eo.message ...
    );
end
try
    ccc.SaveAs([ppath{1} '/' ppath{2} '.ccc']);
catch ne_eo;
    warning( ...
        'neuroelf:ErrorSavingXff', ...
        'Error saving xff of type CCC (%s).', ...
        ne_eo.message ...
    );
end
try
    if tctraw
        ctc.SaveAs([ppath{1} '/' ppath{2} '.ctc']);
    end
catch ne_eo;
    warning( ...
        'neuroelf:ErrorSavingXff', ...
        'Error saving xff of type CTC (%s).', ...
        ne_eo.message ...
    );
end
try
    if tctavg
        act.SaveAs([ppath{1} '/' ppath{2} '.act']);
    end
catch ne_eo;
    warning( ...
        'neuroelf:ErrorSavingXff', ...
        'Error saving xff of type CTC (%s).', ...
        ne_eo.message ...
    );
end
try
    if hasprt
        prt.SaveAs([ppath{1} '/' ppath{2} '.prt']);
    end
catch ne_eo;
    warning( ...
        'neuroelf:ErrorSavingXff', ...
        'Error saving xff of type PRT (%s).', ...
        ne_eo.message ...
    );
end

% clear not needed objects
if nargout < 5
    clearxffobjects(ols((nargout + 1):5));
end
