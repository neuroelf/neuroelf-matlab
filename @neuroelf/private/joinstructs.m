function js = joinstructs(s1, s2, varargin)
% joinstructs  - join structs
%
% FORMAT:       js = joinstructs(s1, s2 [, ...])
%
% Input fields:
%
%       s1, s2      structs
%
% Output fields:
%
%       js          joined struct

% Version:  v0.9c
% Build:    11061315
% Date:     Jun-13 2011, 3:04 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, Jochen Weber
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
if nargin < 2 || ...
   ~isstruct(s1) || ...
   ~isstruct(s2)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end

% linearize
z1 = size(s1);
z2 = size(s2);
s1 = reshape(s1, prod(z1), 1);
s2 = reshape(s2, prod(z2), 1);

% get field names
f1 = fieldnames(s1);
f2 = fieldnames(s2);

% fieldnames the same
if numel(f1) == numel(f2) && ...
    all(strcmp(f1, f2))

    % add together
    js = [s1; s2];

% something needs to be adapted
else

    % multimatch fieldnames
    m1 = multimatch(f1, f2, false, true);
    m2 = multimatch(f2, f1, false, true);
    ad = false;

    % add missing fields
    a1 = 1:numel(f2);
    a1(m1(m1 > 0)) = 0;
    if any(a1 > 0)
        a1 = a1(a1 > 0);
        if isempty(s1)
            s1(1).(f1{1}) = [];
        end
        for fc = 1:numel(a1)
            s1(1).(f2{a1(fc)}) = [];
        end
        if prod(z1) == 0
            s1(:) = [];
        end
        ad = true;
    end
    a2 = 1:numel(f1);
    a2(m2(m2 > 0)) = 0;
    if any(a2 > 0)
        a2 = a2(a2 > 0);
        if isempty(s2)
            s2(1).(f2{1}) = [];
        end
        for fc = 1:numel(a2)
            s2(1).(f1{a2(fc)}) = [];
        end
        if prod(z2) == 0
            s2(:) = [];
        end
        ad = true;
    end

    % fields added?
    if ad
        f1 = fieldnames(s1);
        f2 = fieldnames(s2);
    end

    % different order (or case)
    if ~all(strcmp(f1, f2))

        % non-empty s2
        if ~isempty(s2)

            % resort fieldnames in f2
            m1 = multimatch(f1, f2, false, true);
            c2 = struct2cell(s2);
            s2 = cell2struct(c2(m1, :), f1, 1);
            s2 = reshape(s2, numel(s2), 1);

        % empty s2
        else
            s2 = s1([]);
        end
    end

    % then join
    js = [s1; s2];
end

% reshape
if (numel(s1) < 2 || ...
    findfirst(z1 > 1) > 1) && ...
   (numel(s2) < 2 || ...
    findfirst(z2 > 1) > 1)
    js = reshape(js, ...
        [ones(1, min([findfirst(z1 > 1), findfirst(z2 > 1)]) - 1), numel(js), 1]);
end

% additional arguments
if nargin > 2 && ...
    isstruct(varargin{1})

    % more than one additional struct
    if nargin > 3 && ...
        isstruct(varargin{2})
        js = joinstructs(js, joinstructs(varargin{:}));
    else
        js = joinstructs(js, varargin{1});
    end
end
