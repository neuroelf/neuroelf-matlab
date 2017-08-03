function [minmaxbox] = minmaxbbox(varargin)
% minmaxbbox  - get minimum, maximum bounding box of objects
%
% FORMAT:       minmaxbox = minmaxbbox(v, ...)
%
% Input fields:
%
%       v           either 1x1 VMR or VTC object or cell array
%       ...         further objects
%
% Output fields:
%
%       minmaxbox   4x3 array with minimum (intersect) and
%                   maximum (union) bounding box of objects
%                   valid for .Reframe calls
%
% Note: if all input objects are V16, is returned in natural resolution

% Version:  v0.9a
% Build:    10120112
% Date:     May-17 2010, 10:48 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

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

% argument check
if nargin == 0 || ...
    isempty(varargin{1}) || ...
   ((~iscell(varargin{1}) && ...
     ~isxff(varargin{1}, true)) || ...
    (iscell(varargin{1}) && ...
     ~isxff(varargin{1}{1}, true)))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
tcc = 0;
for ac = 1:nargin
    if ~iscell(varargin{ac})
        varargin{ac} = {varargin{ac}};
    end
    for cc = 1:numel(varargin{ac})
        if numel(varargin{ac}{cc}) ~= 1 || ...
           ~isxff(varargin{ac}{cc}, true)
            error( ...
                'neuroelf:BadArgument', ...
                'Bad argument.' ...
            );
        end
    end
    tcc = tcc + numel(varargin{ac});
end
oc = cell(1, tcc);
tcc = 0;
for ac = 1:nargin
    oc((tcc + 1):(tcc + numel(varargin{ac}))) = varargin{ac}(:);
    tcc = tcc + numel(varargin{ac});
end

% check each type
for cc = tcc:-1:1

    % allowed are AVA/CMP/DDT/GLM/MSK/VDW/VMR/VMP/VTC
    if ~any(strcmpi(oc{cc}.Filetype, ...
         {'ava', 'cmp', 'ddt', 'glm', 'msk', 'vdw', 'vmp', 'vmr', 'vtc'}))
        oc(cc) = [];
    end
end

% reject empty array
if isempty(oc)
    error( ...
        'neuroelf:BadArgument', ...
        'None of the supported types passed in.' ...
    );
end

% start with first
res = zeros(tcc, 3);
fc = oc{1}.BoundingBox.FCube;
minmaxbox = repmat([0, fc(1), fc(1), 0]', [1, 3]);

% iterate from 2 to last
for cc = 1:tcc
    bbox = oc{cc}.BoundingBox;
    nbox = bbox.BBox;
    minmaxbox([1, 4], :) = max(minmaxbox([1, 4], :), nbox);
    minmaxbox([3, 2], :) = min(minmaxbox([3, 2], :), nbox);
    res(cc, :) = bbox.ResXYZ;
end

% check resolution
if ~all(all(res == res(1)))
    error( ...
        'neuroelf:BadArgument', ...
        'Files must match in dimension.' ...
    );
end
