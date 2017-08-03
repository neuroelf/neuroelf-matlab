function vr = dicom_vr(vrtype)
% dicom_vr  - resolve VR values
%
% FORMAT:       vr = dicom_vr(vrtype)
%   - OR -      vrlist = dicom_vr
%
% Input fields:
%
%       vrtype      if given, must be a valid 2-char VR
%
% Output fields:
%
%       vr          1x1 struct with fields
%         .tag      1x2 char, VR tag name
%         .chars    1xN list of valid chars
%         .deblank  false|true (remove leading/trailing spaces/zero bytes)
%         .length   1x2 double describing min/max length
%         .datatype one of 'char', 'int16', 'single', 'uint32', etc.
%
%       vrlist      1x1 struct with VR names as fields and substructs

% Version:  v0.9a
% Build:    11030112
% Date:     May-17 2010, 10:48 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, Jochen Weber
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

% persistent variable
persistent dcm_vr dcm_vr_back;
if isempty(dcm_vr)

    % prepare chars values
    notapplic = '';
    maxchars  = char( 0:255 );
    escchars  = char([27:32:91, 93:255]);
    allchars  = char( 32:126 );
    extchars  = char([32:91, 93:255]);
    defchars  = char([32:91, 93:126]);
    decchars  = char([32, 43, 45:46, 48:57, 69, 101]);
    dtchars   = char([32, 43, 45:46, 48:57]);
    intchars  = char([32, 43, 45, 48:57]);
    codechars = char([32, 48:57, 65:90, 95]);
    datechars = char([46, 48:57]);
    maxlen    = 2^31 - 1;

    % build back
    dcm_vr_back = cell2struct({ ...
        'AE', defchars,               true, [0, 16], 'char'; ...
        'AS', '0123456789DMWYdmwy',  false, [4,  4], 'char'; ...
        'AT', notapplic,             false, [4,  4], 'uint32'; ...
        'CS', codechars,              true, [0, 16], 'char'; ...
        'DA', datechars,             false, [8, 10], 'char'; ...
        'DL', notapplic,             false, [0,  0], 'delim'; ...
        'DS', decchars,               true, [0, 16], 'char'; ...
        'DT', dtchars,                true, [0, 26], 'char'; ...
        'FD', notapplic,             false, [8,  8], 'double'; ...
        'FL', notapplic,             false, [4,  4], 'single'; ...
        'IS', intchars,               true, [0, 12], 'char'; ...
        'LO', allchars,               true, [0, 64], 'char'; ...
        'LT', maxchars,               true, [0, 10240], 'char'; ...
        'OB', notapplic,             false, [0, maxlen], 'uint8'; ...
        'OF', notapplic,             false, [0, maxlen], 'single'; ...
        'OW', notapplic,             false, [0, maxlen], 'uint16'; ...
        'OX', notapplic,             false, [0, maxlen], 'uint816'; ...
        'PN', extchars,               true, [0, 324], 'char'; ...
        'SH', escchars,               true, [0, 16], 'char'; ...
        'SL', notapplic,             false, [4, 4], 'int32'; ...
        'SQ', notapplic,             false, [0, maxlen], 'sequence'; ...
        'SS', notapplic,             false, [0, 2], 'int16'; ...
        'ST', maxchars,               true, [0, 1024], 'char'; ...
        'TM', datechars,              true, [0, 16], 'char'; ...
        'UI', datechars,              true, [0, 64], 'char'; ...
        'UL', notapplic,             false, [4, 4], 'uint32'; ...
        'UN', notapplic,             false, [0, maxlen], 'uint8'; ...
        'US', notapplic,             false, [2, 2], 'uint16'; ...
        'UT', notapplic,             false, [0, maxlen], 'uint8'}, ...
        {'tag', 'chars', 'deblank', 'length', 'datatype'}, 2)';

    % fill structure
    dcm_vr = struct;
    for tc = 1:length(dcm_vr_back)
        dcm_vr.(dcm_vr_back(tc).tag) = dcm_vr_back(tc);
    end
end

% what to return
if nargin < 1
    vr = dcm_vr;
elseif ischar(vrtype) && ...
   ~isempty(vrtype) && ...
    isfield(dcm_vr, upper(vrtype(:)'))
    vr = dcm_vr.(upper(vrtype(:)'));
else
    error( ...
        'xff:BadArgument', ...
        'Invalid DICOM VR argument specified.' ...
    );
end
