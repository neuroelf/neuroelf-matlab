function [xo, pos] = fif_WritePOS(xo, posfile, cspec)
% FIF::WritePos  - write position info for channel coords
%
% FORMAT:       [fif, pos] = fif.WritePos(posfile, cspec)
%
% Input fields:
%
%       posfile     name of POS file
%       cspec       1xC double list of channels to put into POS
%
% Output fields:
%
%       fif         FIF object with loaded headers
%       pos         POS table (Cx9 matrix)

% Version:  v1.1
% Build:    16020310
% Date:     Feb-03 2016, 10:06 AM EST
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

% argument check
if nargin < 3 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'fif') || ...
   ~ischar(posfile) || isempty(posfile) || ~isa(cspec, 'double') || isempty(cspec) || ...
    any(isinf(cspec(:)) | isnan(cspec(:)) | cspec(:) < 1 | cspec(:) ~= fix(cspec(:)))
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
posfile = posfile(:)';
cspec = cspec(:)';

% read in required blocks
fif_ReadInfoHeaders(xo);

% try to get transformation matrix
mat = fif_Value(xo, 'TransformationMatrix');
if isempty(mat)
    error('neuroelf:xff:badFileContent', 'Required transformation matrix not found.');
end
if numel(mat) > 1
    for mc = 1:numel(mat)
        if ischar(mat(mc).From) && ischar(mat(mc).To) && ...
            strcmpi(mat(mc).From, 'device') && strcmpi(mat(mc).To, 'head')
            mat = mat(mc);
            break;
        end
    end
    if numel(mat) > 1
        error('neuroelf:xff:badFileContent', 'Required transformation matrix not found.');
    end
end

% get real matrix
mat = mat.Trans;

% get channel information
chi = fif_Value(xo, 'ChannelInfo');
if isempty(chi)
    error('neuroelf:xff:badFileContent', 'Required channel information not found.');
end
cspec(cspec > numel(chi)) = [];
if isempty(cspec)
    error('neuroelf:xff:badArgument', 'No selected channels in cspec.');
end
NrOfChannels = numel(cspec);

% create pos array
pos = zeros(NrOfChannels, 9);

% iterate over channels
for cc = 1:NrOfChannels
    schi = chi(cspec(cc));
    switch schi.Kind
        case {1, 301}
            if ~iscell(schi.CoilType) || numel(schi.CoilType) < 3
                coilrad = 0.01;
            else
                coilrad = schi.CoilType{3} / 2;
            end
            coilreg = mat * schi.CoilTransformation;
            coildir = coilreg(1:3, 3)';
            coilpt1 = coilreg(1:3, 4)' + coilrad * coilreg(1:3, 1)';
            coilpt2 = coilreg(1:3, 4)' - coilrad * coilreg(1:3, 1)';
            pos(cc, :) = [coilpt1, coilpt2, coildir];
        % case {2}
        otherwise
            error('neuroelf:xff:notYetSupported', ...
                'Getting POS of channels other than MEG not yet supported.');
    end
end

% save POS file
try
    posid = fopen(posfile, 'w');
    if posid < 1
        error('FILE_NOT_WRITABLE');
    end
    for cc = 1:NrOfChannels
        posline = sprintf(' %9.6f  ', pos(cc, :));
        fprintf(posid, '%s\n', posline(1:end-2));
    end
    fclose(posid);
catch xfferror
    error('neuroelf:xff:saveError', 'Error saving POS file: %s.', xfferror.message);
end
