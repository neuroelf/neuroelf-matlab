function ata = transmul(a, b, flag)
% transmul  - return the multiplication of M' * M or M * M2
%
% FORMAT:       Mt_MT_M = transmul(M [, M2 [, flag]])
%
% Input fields:
%
%       M           double matrix (also works along 3rd dim)
%       M2          double matrix
%       flag        determines whether or not to permute first two dims:
%                   - 0: do not permute either M or M2 (default)
%                   - 1: permute dims 1 and 2 in M but not in M2
%                   - 2: permute dims 1 and 2 in M2 but not in M
%                   - 3: permute dims 1 and 2 in both M and M2
%
% Output fields:
%
%       Mt_MT_M     product of transpose with original
%
% Note: this is a MEX (compiled) function!

% Version:  v1.0
% Build:    14091915
% Date:     Sep-19 2014, 3:49 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2014, Jochen Weber
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
   (~isa(a, 'double') && ...
    ~isa(a, 'single'))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing input argument.' ...
    );
end

% return Matlab equivalent if not present
if nargin == 1
    ata = a' * a;

% for additional arguments
else

    % use mtimesnd
    ata = mtimesnd(a, b, flag);
end
