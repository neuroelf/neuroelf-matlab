function ov = partoanalyze75(parfile, pattern)
% partoanalyze75  - convert PAR/REC into Analyze 7.5
%
% FORMAT:       [v =] partoanalyze75(parfile, pattern)
%
% Input fields:
%
%       parfile     PAR filename
%       pattern     output file pattern (e.g. 'RUN1_%03d.img')
%
% Output fields:
%
%       v           SPM vol structure (1xN)
%
% Note: For now, this function uses SPM routines (spm_write_vol)!

% Version:  v0.9b
% Build:    11050712
% Date:     Apr-09 2011, 11:08 PM EST
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

% input check
if nargin < 2 || ...
   ~ischar(parfile) || ...
    isempty(parfile) || ...
    exist(parfile(:)', 'file') ~= 2 || ...
   ~ischar(pattern) || ...
    isempty(pattern) || ...
   ~any(pattern(:)' == '%')
    error( ...
        'SBTools:IllegalUsage', ...
        'Invalid parameters supplied.' ...
    );
end
pattern = pattern(:)';

% read PAR file
parfile = parfile(:)';
try
    pf = readpar(parfile);
    if ~isstruct(pf) || ...
       ~isfield(pf, 'RECData')
        error( ...
            'SBTools:IllegalFile', ...
            'Invalid PAR file name given.' ...
        );
    end
catch ne_eo;
    rethrow(ne_eo);
end

% get data access
rd = pf.RECData;
if isfield(rd, 'Dyn')
    rd = rd.Dyn;
end

% header fields
fov = pf.Parameters.FOV(1:3);
nrs = pf.Parameters.Max_number_of_slices;
rci = find(strcmpi(pf.MatrixHeaders, 'recon_resolution_1'));
rec = pf.MatrixValues(1, [rci, rci + 1]);
sli = find(strcmpi(pf.MatrixHeaders, 'slice_thickness'));
slc = sum(pf.MatrixValues(1, [sli, sli + 1]));
psi = find(strcmpi(pf.MatrixHeaders, 'pixel_spacing_1'));
psp = pf.MatrixValues(1, [psi, psi + 1]);
rfv = [rec, nrs] .* [psp, slc];
if any(abs(rfv - fov) >= 1)
    fov = rfv;
end
dir = fov ./ [rec, nrs];
mat = [ ...
    -dir(1),    0   ,    0   ,  0.5 * (rec(1) + 1) * dir(1); ...
       0   , -dir(2),    0   ,  0.5 * (rec(2) + 1) * dir(2); ...
       0   ,    0   ,  dir(3), -0.5 * (nrs + 1) * dir(3); ...
       0   ,    0   ,    0   ,  1];
v = struct;
v.fname = sprintf(pattern, 1);
v.mat = mat;
v.dim = [rec, nrs];
v.dt = [4, 0];
v.pinfo = [1; 0; 0];
v.n = [1, 1];
v.descrip = pf.Parameters.Protocol_name(1:min(10, numel(pf.Parameters.Protocol_name)));
vy = zeros([rec, nrs]);

% one volume (anat)
if numel(rd) == 1

    for sc = 1:nrs
        vy(:, :, sc) = rd.Slice(sc).IO(:, :);
    end
    ov = spm_write_vol(spm_create_vol(v), vy);

% multi-volume (func)
else

    % iterate over dynamics
    for dc = 1:numel(rd)

        % get slice data
        for sc = 1:nrs
            vy(:, :, sc) = rd(dc).Slice(sc).IO(:, :);
        end

        % write volume
        v.fname = sprintf(pattern, dc);
        if dc == 1
            ov = spm_write_vol(spm_create_vol(v), vy);
            ov = ov(ones(1, numel(rd)));
        else
            ov(dc) = spm_write_vol(spm_create_vol(v), vy);
        end
    end
end
