function xo = fif_LoadChannelData(xo, ctype, cspec)
% FIF::LoadChannelData  - retrieve channel data from FIF file
%
% FORMAT:       fif.LoadChannelData(ctype [, cspec]);
%
% Input fields:
%
%       ctype       Channel type, one of
%                   'MEG', 'RMEG', 'EEG', 'MCG', 'STIM',
%                   'EOG', 'EMG', 'ECG', 'MISC', 'RESP', 'ANY'
%       cspec       Channel specification, if ommited: all channels
%
% No output fields.

% Version:  v1.1
% Build:    16020310
% Date:     Feb-03 2016, 10:11 AM EST
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

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'fif') || ...
   ~ischar(ctype) || ~any(strcmpi(ctype(:)', ...
    {'any', 'ecg', 'eeg', 'emg', 'eog', 'mcg', 'meg', 'misc', 'resp', 'rmeg', 'stim'}))
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end

% shortcut to structure
bc = xo.C;
fst = bc.FIFStructure;

% second argument valid?
if nargin > 1
    % ...
else

    % find and read channel info elements
    cspec = find(fst.Lookup == 203);

    % read those
    fif_ReadElements(xo, cspec);
end
