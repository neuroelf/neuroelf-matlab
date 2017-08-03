function l = ti_renderlayer(l, h, w)
% transimg::PRIVATE::ti_renderlayer  - render a non-full layer to buffer
%
% FORMAT:       l = ti_renderlayer(l, h, w)
%
% Input fields:
%
%       l           layer struct (see transimg class definition)
%       h, w        height and width
%
% Output fields:
%
%       l           struct with updated RPixel, RAlpha, and IsRendered
%
% This is an interally called function!
%
% Using flexinterpn_method, spmtrf.

% Version:  v1.0
% Build:    16021216
% Date:     Feb-12 2016, 4:39 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2012 - 2016, Jochen Weber
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

% neuroelf library
global ne_methods;

% argument check
if nargin ~= 3 || ~isa(h, 'double') || numel(h) ~= 1 || isinf(h) || isnan(h) || h <= 0 || h ~= fix(h) || ...
   ~isa(w, 'double') || numel(w) ~= 1 || isinf(w) || isnan(w) || w <= 0 || w ~= fix(w) || ...
   ~isstruct(l) || numel(l) ~= 1 || numel(fieldnames(l)) ~= 8 || ...
   ~all(strcmp(fieldnames(l), {'Type'; 'Pixel'; 'Alpha'; 'Trans'; 'IsRendered'; 'RPixel'; 'RAlpha'; 'Ref'})) || ...
   ~ischar(l.Type) || numel(l.Type) ~= 1 || ~any(l.Type == 'fpstx')
    error('neuroelf:transimg:badArgument', 'Invalid layer structure.');
end

% nothing to do
if any(l.Type == 'fp') || l.IsRendered
    return;
end

% what to do
switch (l.Type)

    % shape
    case 's'

        % to be implemented
        error('neuroelf:transimg:notYetImplemented', 'Shapes not yet implemented.');

    % transformation
    case 't'

        % depending on transformation type
        t = l.Trans;

        % translation and rotation
        if numel(t) == 3
            t = ne_methods.spmtrf([0, t(1:2)], [t(3), 0, 0]);

        % small matrix
        elseif isequal(size(t), [3, 3])
            t = [[1, 0, 0, 0]; [zeros(3, 1), t]];

        % unsupported
        elseif ~isequal(size(t), [4, 4])
            error('neuroelf:transimg:badTransformation', 'Invalid transformation information.');
        end

        % get flexinterpn argument
        fic = [inf, inf, inf; 1, 1, 1; 1, 1, 1; size(l.Pixel, 3), h, w];

        % apply transformation to pixel
        if size(l.Pixel, 3) > 1
            l.RPixel = permute(ne_methods.flexinterpn_method(...
                permute(l.Pixel, [3, 1, 2]), fic, t, 'cubic'), [2, 3, 1]);
        else
            l.RPixel = reshape(ne_methods.flexinterpn_method(...
                reshape(l.Pixel, [1, size(l.Pixel)]), fic, t, 'cubic'), [h, w]);
        end

        % alpha?
        if numel(l.Alpha) ~= 1
            fic(3, 1) = 1;
            l.RAlpha = reshape(ne_methods.flexinterpn_method(...
                reshape(l.Alpha, [1, size(l.Alpha)]), fic, t, 'cubic'), [h, w]);
        else
            l.RAlpha = l.Alpha;
        end

        % rendered
        l.IsRendered = true;

    % text
    case 'x'

        % to be implemented
        error('neuroelf:transimg:notYetImplemented', ...
            'Text layers not yet implemented.');
end
