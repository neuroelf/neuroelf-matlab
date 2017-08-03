function mdm = glm_GenerateMDM(xo, prtlup)
% GLM::GenerateMDM  - return MDM with names of files and settings
%
% FORMAT:       mdm = glm.GenerateMDM([prtlup]);
%
% Input fields:
%
%       prtlup      either of {'none'}, 'sdm', or 'xtc'
%
% Output fields:
%
%       mdm         MDM object with filenames and settings

% Version:  v1.1
% Build:    16021613
% Date:     Feb-16 2016, 1:56 PM EST
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

% check arguments
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'glm')
    error('neuroelf:xff:badArgument', 'Invalid object handle in call.');
end
if nargin < 2 || ~ischar(prtlup) || isempty(prtlup) || ~any(lower(prtlup(1)) == 'nsx')
    prtlup = 'n';
else
    prtlup = lower(prtlup(1));
end
bc = xo.C;
s = bc.Study(:);

% generate mdm
mdm = xff('new:mdm');

% get content
mdmc = mdm.C;

% make global settings
switch (bc.ProjectType)
    case 0
        mdmc.TypeOfFunctionalData = 'FMR';
    case 1
        mdmc.TypeOfFunctionalData = 'VTC';
    case 2
        mdmc.TypeOfFunctionalData = 'MTC';
end
mdmc.RFX_GLM = double(bc.ProjectTypeRFX > 0);
mdmc.PSCTransformation = double(bc.TransformationType == 3);
mdmc.zTransformation = double(bc.TransformationType == 1);
mdmc.SeparatePredictors = bc.SeparatePredictors;
mdmc.NrOfStudies = numel(s);

% generate list of files
if bc.ProjectType <= 1
    xtc_rtc = cell(numel(s), 2);
else
    xtc_rtc = cell(numel(s), 3);
end
for sc = 1:numel(s)
    if bc.ProjectType > 1
        if ischar(s(sc).NameOfSSMFile)
            xtc_rtc{sc, 1} = s(sc).NameOfSSMFile(:)';
        else
            xtc_rtc{sc, 1} = '';
        end
    end
    if ischar(s(sc).NameOfAnalyzedFile)
        xtc_rtc{sc, end-1} = s(sc).NameOfAnalyzedFile(:)';
    else
        xtc_rtc{sc, end-1} = '';
    end
    if ischar(s(sc).NameOfSDMFile)
        xtc_rtc{sc, end} = s(sc).NameOfSDMFile(:)';
    else
        xtc_rtc{sc, end} = '';
    end
end

% PRT lookup
switch (prtlup)

    % SDM name based lookup
    case 's'

        % loop over studies
        for sc = 1:numel(s)

            % if PRT with same name exists
            if numel(xtc_rtc{sc, end}) > 4 && ...
                exist([xtc_rtc{sc, end}(1:end-4), '.prt'], 'file') == 2
                xtc_rtc{sc, end} = [xtc_rtc{sc, end}(1:end-4), '.prt'];
            end
        end

    % time-course object based lookup
    case 'x'

        % loop with error handling
        for sc = 1:numel(s)
            try
                soc = xff(xtc_rtc{sc, end-1}, 'h');
                soc = soc.C;
                if iscell(soc.NameOfLinkedPRT) && ~isempty(soc.NameOfLinkedPRT)
                    xtc_rtc{sc, end} = so.NameOfLinkedPRT{1}(:)';
                elseif ~isempty(soc.NameOfLinkedPRT)
                    xtc_rtc{sc, end} = so.NameOfLinkedPRT(:)';
                end
            catch xfferror
                neuroelf_lasterr(xfferror);
            end
        end
end

% set in mdm
mdmc.XTC_RTC = xtc_rtc;
mdm.C = mdmc;
