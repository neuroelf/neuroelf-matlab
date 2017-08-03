function h = scaleimage(cdata, varargin)
% scaleimage  - automatically scale image data for image
%
% FORMAT:       [h = ] scaleimage(C, ...)
%
% Input fields:
%
%       C           MxN image data that is shown with a 256 grayscale
%       ...         additional arguments passed to image
%
% Output fields:
%
%       h           image handle

% Version:  v1.0
% Build:    15040310
% Date:     Apr-03 2015, 10:01 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2015, Jochen Weber
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

% persistent last handle
persistent i_scilasth;
if isempty(i_scilasth)
    i_scilasth = -1;
end

% argument check
if ~isnumeric(cdata) || ...
    ndims(cdata) > 2
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid color data given.' ...
    );
end
mmm = minmaxmean(cdata, 4);
cds = floor(1 + (255.99 / (eps + mmm(2) - mmm(1))) .* (double(cdata) - mmm(1)));
if nargin < 3 || ...
   ~ischar(varargin{end-1}) || ...
   ~any(strcmpi(varargin{end-1}(:)', {'image', 'parent'})) || ...
   (~isa(varargin{end}, 'double') && ...
    ~isa(varargin{end}, 'matlab.graphics.axis.Axes') && ...
    ~isa(varargin{end}, 'matlab.graphics.primitive.Image')) || ...
    numel(varargin{end}) ~= 1 || ...
   ~ishandle(varargin{end}) || ...
   ~any(strcmpi(get(varargin{end}, 'Type'), {'axes', 'image'}))
    if ishandle(i_scilasth) && ...
       ~isempty(get(0, 'CurrentFigure')) && ...
       ~isempty(get(get(0, 'CurrentFigure'), 'CurrentAxes')) && ...
        any(get(get(get(0, 'CurrentFigure'), 'CurrentAxes'), 'Children') == i_scilasth)
        h = i_scilasth;
        set(h, 'CData', cds);
    else
        if ~isempty(get(0, 'CurrentFigure')) && ...
            isempty(get(get(0, 'CurrentFigure'), 'Children'))
            h = get(0, 'CurrentFigure');
        else
            h = figure;
        end
        h = axes('Parent', h);
        h = image(cds, varargin{:}, 'Parent', h);
        i_scilasth = h;
        colormap((0:1/255:1)' * ones(1,3));
    end
else
    if strcmpi(get(varargin{end}, 'Type'), 'image')
        h = varargin{end};
        set(h, 'CData', cds);
    else
        h = image(cds, varargin{:});
        i_scilasth = h;
    end
end
drawnow;
