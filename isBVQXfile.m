function tf = isBVQXfile(varargin)
%ISBVQXFILE  Perform object verification (for compatibility with old code).
%   TF = ISBVQXFILE(OBJ) returns true if the input is of type XFF.
%
%   TF = ISBVQXFILE(OBJ, true) returns false if the input type (class) is
%   not XFF, and a sized true/false array for each array member depending
%   on whether OBJ(INDEX) is valid.
%
%   TF = ISBVQXFILE(OBJ, TYPE) returns false if the input type (class) is
%   not XFF, and a sized true/false array for each array member depending
%   on whether OBJ(INDEX) is a valid XFF object of type TYPE.
%
%   See also XFF

% Version:  v1.1
% Build:    16031110
% Date:     Mar-11 2016, 10:05 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2016, Jochen Weber
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

% base argument check
if nargin < 1 || nargin > 2
    error('neuroelf:xff:badNumberOfInputs', 'Invalid number of inputs.');
end

% return false (otherwise, it uses the class method!)
tf = false;
