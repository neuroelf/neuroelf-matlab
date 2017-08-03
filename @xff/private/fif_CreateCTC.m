function [ctc, stimchan] = fif_CreateCTC(xo, ctcfile, cspec)
% FIF::CreateCTC  - create CTC object from FIF data
%
% FORMAT:       [ctc, stimchan] = fif.CreateCTC(ctcfile, [cspec])
%
% Input fields:
%
%       ctcfile     name of CTC output file
%       cspec       1x1 struct with optional fields:
%        .CoilType  list of coiltypes from channel info to export
%                   (see fif.Value('ChannelInfo') for details)
%        .ExpList   1xN list of channels to export, e.g. [1:3:318]
%        .StimAuto  if given, try to guess stimulation channels
%        .StimList  1xN list of stimulus channels (to import into PRT)
%
% Output fields:
%
%       ctc         the created object (with transio CTCData!)
%       stimchan    data of stimulus channels (if specified with nonstim)
%
% Note: if the cspec is not given, by default all channels are exported
%
% Using: fifio.

% Version:  v1.1
% Build:    16020310
% Date:     Feb-03 2016, 10:22 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2011, 2014, 2016, Jochen Weber
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

% neuroelf library and factory
global ne_methods xffsngl;

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'fif') || ...
   ~ischar(ctcfile) || isempty(ctcfile)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
ctcfile = ctcfile(:)';
if nargin < 3 || ~isstruct(cspec) || numel(cspec) ~= 1
    cspec = struct;
end
if isfield(cspec, 'CoilType')
    if ~isa(cspec.CoilType, 'double') || ...
        max(size(cspec.CoilType)) ~= numel(cspec.CoilType) || ...
        any(isinf(cspec.CoilType) | isnan(cspec.CoilType) | ...
            cspec.CoilType < 1 | cspec.CoilType ~= fix(cspec.CoilType))
        cspec = rmfield(cspec, 'CoilType');
    end
end
if isfield(cspec, 'ExpList')
    if ~isa(cspec.ExpList, 'double') || ...
        max(size(cspec.ExpList)) ~= numel(cspec.ExpList) || ...
        any(isinf(cspec.ExpList) | isnan(cspec.ExpList) | ...
            cspec.ExpList < 1 | cspec.ExpList ~= fix(cspec.ExpList))
        cspec = rmfield(cspec, 'ExpList');
    end
end
if isfield(cspec, 'StimList')
    if ~isa(cspec.Stimlist, 'double') || ...
        max(size(cspec.StimList)) ~= numel(cspec.StimList) || ...
        any(isinf(cspec.StimList) | isnan(cspec.StimList) | ...
            cspec.StimList < 1 | cspec.StimList ~= fix(cspec.StimList))
        cspec = rmfield(cspec, 'StimList');
    end
end

% make sure the info blocks are read
fif_ReadInfoHeaders(xo);

% get content
bc = xo.C;

% get fif shortcut and check for meta information
fif = bc.FIFStructure;
if ~any(fif.BlockLookup(1, :) == 101)
    error('neuroelf:xff:badInputFile', 'FIF file is missing the MetaInformation block.');
end

% try to read NrOfChannels and SamplingFrequency
NrOfChannels = -1;
SamplingFrequency = -1;
nci = find(fif.Lookup == 200);
sfi = find(fif.Lookup == 201);
if isempty(nci) || isempty(sfi)
    error('neuroelf:xff:badInputFile', 'No NrOfChannels and/or SamplingFrequency in FIF.');
end

% lookup blocks of tags
fifio = ne_methods.fifio;
for ncic = nci
    eblock = fifio(fif, 'elemblock', ncic);
    if ~isempty(eblock) && eblock(1) == 101
        NrOfChannels = fif.Children(ncic).Value;
        break;
    end
end
for sfic = sfi
    eblock = fifio(fif, 'elemblock', sfic);
    if ~isempty(eblock) && eblock(1) == 101
        SamplingFrequency = fif.Children(sfic).Value;
        break;
    end
end
if NrOfChannels < 0 || SamplingFrequency < 0
    error('neuroelf:xff:badInputFile', ...
        'Unable to determine NrOfChannels and/or SamplingFrequency.');
end

% get channel info
chinfo = fif_Value(xo, 'ChannelInfo');
if numel(chinfo) ~= NrOfChannels
    error('neuroelf:xff:badInputFile', ...
        'NrOfChannels from FIF header does not match number of Infos.');
end

% find out about layout in FIF file
dbi = find(fif.Lookup == 300);
blockm = find(diff(dbi) - 1) + 1;
if ~isempty(blockm)
    error('neuroelf:xff:badInputFile', ...
        'Multiple blocks of data buffers in FIF file unsupported.');
end
buffnum = numel(dbi);

% check first and last buffer for size
fif = fifio(fif, 'readelem', [dbi(1), dbi(end)]);
fbsize = numel(fif.Children(dbi(1)).Value);
lbsize = numel(fif.Children(dbi(end)).Value);
if (fbsize / NrOfChannels) ~= fix(fbsize / NrOfChannels) || ...
   (lbsize / NrOfChannels) ~= fix(lbsize / NrOfChannels)
    error('neuroelf:xff:badInputFile', ...
        'Size of data buffer elements and NrOfChannels don''t match.');
end
sperbuff = fbsize / NrOfChannels;

% stimulus channel specification
exportc = 1:NrOfChannels;
stimulc = [];
if isfield(cspec, 'StimAuto')
    stich = fif_Value(xo, 'StimulusChannels');
    if ~isempty(stich)
        if iscell(stich)
            for sc = 1:numel(stich)
                stimulc = union(stimulc, stich{sc});
            end
        else
            stimulc = stich(:)';
        end
    end
    if isfield(cspec, 'StimList')
        cspec = rmfield(cspec, 'StimList');
    end
end
if isfield(cspec, 'StimList')
    stimulc = cspec.StimList(:)';
end
if ~isempty(stimulc)
    exportc = setdiff(exportc, stimulc);
end
stimulc(stimulc > NrOfChannels) = [];

% export channel specification
if isfield(cspec, 'CoilType')
    ccoiltype = {chinfo.CoilType};
    coiltype = false(1, NrOfChannels);
    for cc = 1:NrOfChannels
        if any(cspec.CoilType == ccoiltype{cc})
            coiltype = true;
        end
    end
    exportc = find(coiltype);
elseif isfield(cspec, 'ExpList')
    exportc = cspec.ExpList(:)';
end
exportc(exportc > NrOfChannels) = [];

% get numbers
NrOfChannelsEx = numel(exportc);
NrOfChannelsSt = numel(stimulc);
if NrOfChannelsSt > 0
    stimchan = single(0);
    stimchan(buffnum * fbsize * NrOfChannelsSt / NrOfChannels) = 0;
    stimchan = reshape(stimchan, [NrOfChannelsSt, numel(stimchan) / NrOfChannelsSt]);
else
    stimchan = single(zeros(0, buffnum * fbsize / NrOfChannels));
end
NrOfSamples = (buffnum - 1) * sperbuff + lbsize / NrOfChannels;

% create CTC in memory
ctc = xff('new:ctc');
ctco = ctc.C;

% set values in CTC and write bogus
ctco.SampleOrdering = 2;
ctco.NrOfChannels = NrOfChannelsEx;
ctco.NrOfSamples = 1024;
ctco.SamplingFrequency = double(SamplingFrequency);
ctco.CTCData = single(zeros(NrOfChannelsEx, 1024));
ctc.C = ctco;
try
    aft_SaveAs(ctc, ctcfile);
catch xfferror
    delete(ctc);
    error('neuroelf:xff:saveFailed', ...
        'Error saving CTC file to ''%s'': %s.', ctcfile, xfferror.message);
end
delete(ctc);

% get and re-set transiosize for ctc
ctcidx = xffsngl.EXT.ctc{2};
ctctiosize = xffsngl.BFF(ctcidx).TransIOSize;
try
    xffsngl.BFF(ctcidx).TransIOSize = 256;
    xffsngl.FF.bff(ctcidx).TransIOSize = 256;
    ctc = xff(ctcfile);
    ctco = ctc.C;
    tiostr = ctco.CTCData;
    xffsngl.BFF(ctcidx).TransIOSize = ctctiosize;
    xffsngl.FF.bff(ctcidx).TransIOSize = ctctiosize;
catch xfferror
    xffsngl.BFF(ctcidx).TransIOSize = ctctiosize;
    xffsngl.FF.bff(ctcidx).TransIOSize = ctctiosize;
    delete(ctc);
    error('neuroelf:xff:openFailed', ...
        'Error opening CTC in transio mode: %s.', xfferror.message);
end

% now set the true values and expand file!
ctco.NrOfSamples = NrOfSamples;

% expand file
try
    tio = transio(tiostr.FileName, 'ieee-le', 'single', tiostr.IOOffset, ...
        [NrOfChannelsEx, NrOfSamples], 1);

    % and also set NrOfSamples (a bit brute force, but we know the header!)
    tio2 = transio(tiostr.FileName, 'ieee-le', 'uint32', 20, [1, 1]);
    tio2(1) = uint32(NrOfSamples);
    if double(tio2(1)) ~= NrOfSamples
        error('TRANSIO_FAILURE');
    end
catch xfferror
    delete(ctc);
    error('neuroelf:xff:transIOError', ...
        'Error expanding file for CTC data: %s.', xfferror.message);
end
ctco.CTCData = tio;
ctc.C = ctco;

% fill transio object
spc = 1;
for buc = 1:buffnum

    % read buffer
    fif = fifio(fif, 'readelem', dbi(buc));

    % put data into tio
    buffdata = fif.Children(dbi(buc)).Value;
    buffdata = reshape(buffdata, [NrOfChannels, numel(buffdata) / NrOfChannels]);
    nrofsamp = size(buffdata, 2);
    tio(:, spc:(spc + nrofsamp - 1)) = buffdata(exportc, :);

    % keep stimulus data ?
    if NrOfChannelsSt > 0
        stimchan(:, spc:(spc + nrofsamp - 1)) = buffdata(stimulc, :);
    end
    spc = spc + nrofsamp;

    % free memory again
    fif = fifio(fif, 'freeelem', dbi(buc));
end
