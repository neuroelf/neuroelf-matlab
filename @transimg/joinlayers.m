function ti = joinlayers(ti, l, varargin)
% transimg::joinlayers  - join layers
%
% Using: joinlayers.

% Version:  v0.9d
% Build:    14082618
% Date:     Aug-26 2014, 6:03 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2012, 2014, Jochen Weber
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

% global variables for storage
global tiobj ...
       tiobjlup ...
       ne_methods;

% check arguments
lup = find(tiobjlup == ti.L);
if numel(lup) ~= 1
    error( ...
        'transimg:ObjectRemoved', ...
        'Object removed from global storage.' ...
    );
end
if nargin < 2 || ...
   ~isa(l, 'double') || ...
    any(isinf(l(:)) | isnan(l(:)) | l(:) < 1)
    error( ...
        'transimg:BadArgument', ...
        'Invalid layer selection.' ...
    );
end
l = floor(l(:));
l(l > numel(tiobj(lup).Layer)) = [];
if isempty(l)
    error( ...
        'transimg:BadArgument', ...
        'Invalid layer selection.' ...
    );
end
if numel(unique(l)) == 1
    return;
end
l1 = l(1);

% types must all be full or partial
lt = {tiobj(lup).Layer(l).Type};
if any(strcmp(lt, 'p'))
    isp = true;
    if ~all(strcmp(lt, 'p'))
        error( ...
            'transimg:InvalidCombination', ...
            'All joined layers must be either partial or full.' ...
        );
    end
    lsz = zeros(numel(l), 4);
    for lc = 1:numel(lt)
        tlay = tiobj(lup).Layer(l(lc));
        lsz(lc, :) = [size(tlay.Pixel, 1), size(tlay.Pixel, 2), tlay.Trans(1:2)];
    end
    if any(any(diff(lsz, 1, 1)))
        error( ...
            'transimg:InvalidCombination', ...
            'Partial layers must match in offset and size to be joined.' ...
        );
    end
else
    isp = false;
end

% render individual layers
lr = {tiobj(lup).Layer(l).IsRendered};
for lc = 1:numel(lt)
    if ~ischar(lt{lc}) || ...
       ~any(strcmp(lt{lc}, {'f', 'p', 's', 't', 'x'}))
        error( ...
            'transimg:BadLayerType', ...
            'Invalid layer type in layer %d: %s.', ...
            l(lc), lt{lc} ...
        );
    end
    if any(lt{lc} == 'fp') || ...
        lr{lc}
        continue;
    end
    tiobj(lup).Layer(l(lc)) = ti_renderlayer(tiobj(lup).Layer(l(lc)), tiobj(lup).Height, tiobj(lup).Width);
end

% put back into .Rendered
if isp
    ttiobj = tiobj(lup);
    ttiobj.Height = size(ttiobj.Layer(l(1)).Pixel, 1);
    ttiobj.Width = size(ttiobj.Layer(l(1)).Pixel, 2);
    [p, a] = ne_methods.joinlayersc(ttiobj, l, varargin{:});
else
    [p, a] = ne_methods.joinlayersc(tiobj(lup), l, varargin{:});
end
tiobj(lup).Layer(l1).Pixel = p;
tiobj(lup).Layer(l1).Alpha = a;
l = unique(l);
l(l == l1) = [];
tiobj(lup).Layer(l) = [];
tiobj(lup).IsRendered = false;
