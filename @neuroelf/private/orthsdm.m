function orthsdm(filename, newfilename)
% orthsdm  - orthogonalize SDM object in file
%
% FORMAT:       orthsdm(filename [, newfilename]);
%
% Input fields:
%
%       filename    SDM file
%       newfilename SaveAs filename (default: *_orth.sdm)
%
% No output fields

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

% argument check
if nargin < 1 || ...
   ~ischar(filename) || ...
    isempty(filename)
    error( ...
        'neuroelf:BadArgument', ...
        'Valid matrix of vectors without Inf/Nans are required.' ...
    );
end
filename = filename(:)';
if nargin < 2 || ...
   ~ischar(newfilename) || ...
    isempty(newfilename)
    [fp{1:3}] = fileparts(filename);
    if isempty(fp{1})
        fp{1} = '.';
    end
    newfilename = [fp{1} filesep fp{2} '_orth.sdm'];
end

% try loading SDM
sdm = [];
sdml = true;
try
    [sdm, sdml] = xff(filename);
    if ~isxff(sdm, 'sdm')
        error('NO_SDM');
    end
catch ne_eo;
    if isxff(sdm, true) && ...
        sdml
        sdm.ClearObject;
    end
    error( ...
        'neuroelf:BadArgument', ...
        'No SDM file loaded: %s.', ...
        ne_eo.message ...
    );
end

% orthogonalize SDMMatrix
if sdm.FirstConfoundPredictor < 1 || ...
    sdm.FirstConfoundPredictor > size(sdm.SDMMatrix, 2)
    sdm.SDMMatrix = orthvecs(sdm.SDMMatrix);
else
    sdm.SDMMatrix(:, 1:sdm.FirstConfoundPredictor-1) = ...
        orthvecs(sdm.SDMMatrix(:, 1:sdm.FirstConfoundPredictor-1));
end

% SaveAs
sdm.SaveAs(newfilename);
if sdml
    sdm.ClearObject;
end
