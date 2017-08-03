function [ctc, prt, stimdata] = importfif(outfile, cspec, fiffile)
% importfif  - import FIF file
%
% FORMAT:       [ctc, prt, sdata] = importfif(outfile, cspec [, fiffile])
%
% Input fields:
%
%       outfile     CTC output filename (may include path)
%       cspec       condition specification (used for
%                   PRT::ImportStimChannels, if not given, no PRT created)
%       fiffile     FIF filename, interactively if not given or empty
%
% Output fields:
%
%        ctc        CTC object
%        prt        PRT object
%        sdata      stimulation data from FIF::CreateCTC method

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:20 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2016, Jochen Weber
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
   ~ischar(outfile) || ...
    isempty(outfile)
    error( ...
        'neuroelf:BadArguments', ...
        'Bad or missing argument.' ...
    );
end
outfile = outfile(:)';
[outp{1:3}] = fileparts(outfile);
if isempty(outp{1})
    outp{1} = '.';
end
if ~strcmpi(outp{3}, '.ctc')
    outfile = [outp{1} '/' outp{2} '.ctc'];
end
if nargin < 2 || ...
   ~isstruct(cspec) || ...
    isempty(cspec) || ...
   ~isfield(cspec, 'Name') || ...
   ~isfield(cspec, 'Color') || ...
   ~isfield(cspec, 'Pattern') || ...
   ~isfield(cspec, 'Duration')
    cspec = [];
else
    cspec = cspec(:)';
end
if nargin < 3 || ...
   ~ischar(fiffile) || ...
    isempty(fiffile) || ...
    exist(fiffile(:)', 'file') ~= 2
    fiffile = '*.fif';
end

% set transiosize of CTC
fif = xff(fiffile);
fif.ReadInfoHeaders;

% get number of channels
noc = fif.Value('NrOfSensorChannels');

% fill channel specification
chspec = struct;
chspec.ExpList = sort([1:3:noc, 2:3:noc]);
chspec.StimAuto = true;

% create CTC and get stim data
[ctc, stimdata] = fif.CreateCTC(outfile, chspec);

% write positional data
fif.WritePOS([outp{1} '/' outp{2} '.pos'], chspec.ExpList);
fif.WriteSFH([outp{1} '/' outp{2} '.sfh']);

% FIF no longer needed
fif.ClearObject;

% create protocol
if ~isempty(cspec)
    prt = xff('new:prt');

    % import stimulus channels
    prt.ImportStimChannels(stimdata, ctc.SamplingFrequency, cspec);

    % save PRT
    prt.SaveAs([outp{1} '/' outp{2} '.prt']);

    % keep protocol ?
    if nargout < 2
        prt.ClearObject;
    end

% no protocol
else
    prt = [];
end

% keep CTC open
if nargout < 1
    ctc.ClearObject;
end
