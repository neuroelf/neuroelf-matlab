function ti = transimg(varargin)
% transimg  - class and constructor for transimg objects
%
% Fields (properties):
%
%  .Width                           width (1x1 double)
%  .Height                          height (1x1 double)
%  .Background                      background RGB (1x3 uint8)
%  .IsRendered                      boolean flag, needs rendering?
%  .Rendered                        after call to render contains the
%                                   Width x Height x 3 RGB uint8 image
%  .Layer(L).Pixel                  pixel information of layer L
%           .Alpha                  alpha information of layer L
%  .Handle                          associated image handle (or empty)
%
% Methods
%
% TI = transimg(WIDTH, HEIGHT); constructor
% addlayer(TI, PIXELS, ALPHA);  add a layer of pixels/alpha information
% avglayers(TI, AVGINDEX);      average layers (into first of index)
% delete(TI);                   remove object from global storage
% dellayer(TI, DELINDEX);       remove the named layer(s)
% display(TI);                  set handle with CData (create if necessary)
% joinlayers(TI, JOININDEX);    joins the named layers into one layer
% render(TI);                   render the layers into .Rendered
% sethandle(TI, HANDLE);        set graphics handle (of type image)
%
% Using: emptystruct, findfirst.

% Version:  v1.1
% Build:    16052021
% Date:     May-20 2016, 9:14 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2012, 2014, 2016, Jochen Weber
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
global tiobj tiobjlup ne_methods;
   
% initialize global variable (only once)
if isempty(tiobj)

    % neuroelf library required (init)
    nelf = neuroelf;

    % then generate image
    tiobj = newtransimg;
    tiobjlup = 0;
end

% if input is of type transimg, special uses
if nargin > 0 && isa(varargin{1}, 'transimg')

    % WHICH ARE NOT YET IMPLEMENTED
end

% preset width and height to 0
height = 0;
width = 0;

% lookup constructor
if nargin == 1 && isa(varargin{1}, 'double') && numel(varargin{1}) == 1
    lup = ne_methods.findfirst(tiobjlup == varargin{1});
    if ~isempty(lup)
        ti = class(struct('L', varargin{1}), 'transimg');
        return;
    end
    error( ...
        'transimg:BadLookup', ...
        'Lookup of internal handle failed.' ...
    );
end

% image constructor
if nargin > 1 && isa(varargin{1}, 'uint8') && numel(varargin{1}) > 1 && ...
   (isa(varargin{2}, 'double') || isa(varargin{2}, 'single')) && ...
    size(varargin{1}, 1) == size(varargin{2}, 1) && ...
    size(varargin{1}, 2) == size(varargin{2}, 2) && ...
    ndims(varargin{2}) == 2
    try
        ti = transimg(size(varargin{1}, 2), size(varargin{1}, 1), varargin{3:end});
        addlayer(ti, varargin{1}, varargin{2});
        return;
    catch ne_eo;
        rethrow(ne_eo);
    end
end

% width and height constructor
if nargin > 1 && isa(varargin{1}, 'double') && numel(varargin{1}) == 1 && ...
    isa(varargin{2}, 'double') && numel(varargin{2}) == 1
    width = varargin{1};
    height = varargin{2};
end
if nargin > 2 && (isa(varargin{3}, 'double') || isa(varargin{3}, 'single') || ...
     isa(varargin{3}, 'uint8')) && numel(varargin{3}) == 3
    bgcol = double(varargin{3}(:)');
    bgcol(isinf(bgcol) | isnan(bgcol) | bgcol < 0) = 0;
    while any(bgcol > 1)
        bgcol = (1 / 255) .* bgcol;
    end
else
    bgcol = [0, 0, 0];
end
bgcol = uint8(round(255 .* bgcol));
uo = uint8(ones(height, width));

% create new image im global storage
tiobj(end + 1) = newtransimg;
lup = rand(1, 1);
while any(tiobjlup == lup)
    lup = rand(1, 1);
end
tiobjlup(end + 1) = lup;

% make settings
tiobj(end).Width = width;
tiobj(end).Height = height;
tiobj(end).Background = bgcol;
tiobj(end).Rendered = cat(3, bgcol(1) .* uo, bgcol(2) .* uo, bgcol(3) .* uo);
tiobj(end).IsRendered = true;
tiobj(end).Handle = [];

% then create object
ti = class(struct('L', lup), 'transimg');

% subfunction
function no = newtransimg
global ne_methods;
no = struct( ...
    'Width'     , 0, ...
    'Height'    , 0, ...
    'Background', [0, 0, 0], ...
    'IsRendered', false, ...
    'Rendered'  , uint8(zeros(0, 0, 3)), ...
    'Layer'     , ne_methods.emptystruct({'Type', 'Pixel', 'Alpha', 'Trans', 'IsRendered', 'RPixel', 'RAlpha', 'Ref'}), ...
    'Handle'    , []);
