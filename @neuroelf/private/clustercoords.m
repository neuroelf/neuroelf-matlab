function [cc, cl, cs, cv] = clustercoords(vol, meth, kthr)
% clustercoords  - get list of clusters in volume
%
% FORMAT:       [cc, cl, cs, cv] = clustercoords(vol [, meth [, kthr]])
%
% Input fields:
%
%       vol         XxYxZ logical array, true for indices to be clustered
%       meth        either of 'face', 'edge', or 'vertex' (default: face)
%       kthr        cluster-size threshold (default: 1)
%
% Output fields:
%
%       cc          cluster coords, Cx4 list of coordinates, where
%                   the 4th column is a cluster index within 1...N
%       cl          1xN cell array with Px3 coordinates
%       cs          sizes of clusters
%       cv          clustered volume (see clustercoordsc)
%
% Note: this now uses the compiled function clustercoordsc !!!

% Version:  v0.9b
% Build:    11051611
% Date:     Apr-09 2011, 1:55 PM EST
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
   ~islogical(vol) || ...
    isempty(vol)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin < 2 || ...
   ~ischar(meth) || ...
    isempty(meth)
    meth = 'f';
else
    meth = lower(meth(1));
end
if ~any('efv' == meth)
    meth = 'f';
end
if nargin < 3 || ...
   ~isa(kthr, 'double') || ...
    numel(kthr) ~= 1 || ...
    isinf(kthr) || ...
    isnan(kthr) || ...
    kthr < 1 || ...
    kthr > numel(vol)
    kthr = 1;
else
    kthr = round(kthr);
end
switch (meth)
    case {'e'}
        meth = 2;
    case {'f'}
        meth = 1;
    otherwise
        meth = 3;
end

% check for compiled function
try
    [cs, cv, cc, cl] = clustercoordsc(vol, meth, kthr);
catch ne_eo;
    rethrow(ne_eo);
end
