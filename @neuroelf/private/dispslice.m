function [varargout] = dispslice(img, h, varargin)
% dispslice  - displays a 2-D dataset as an image
%
% FORMAT:       dispslice(slicedata [,handle])
%     or        dispslice(figcaption [,slicedim])
%
% Input fields:
%
%       slicedata   NxN double image data (2-D)
%       handle      valid image handle to use
%       slicedim    dimension for display

% Version:  v0.9d
% Build:    14061315
% Date:     Jun-13 2014, 3:00 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2014, Jochen Weber
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

% persistent variable
persistent i_myfig;

% enough arguments ?
if nargin < 1
    error( ...
        'neuroelf:TooFewArguments',...
        'Too few arguments. Try ''help %s''.',...
        mfilename ...
    );
end

% what kind of input
if ischar(img)

    % parse other arguments
    try

        % if no more input is given, set handle to empty
        if nargin == 1
            h = [];
        end

        % is handle is a char number, try conversion
        if ischar(h)
            h = str2double(h);
        end

        % what is the first argument
        switch(img),

            % axes access
            case {'axes'}
                varargout{1} = i_myfig.a;
                return;

            % closefig request
            case {'closefig'}
                delete(i_myfig.f);
                i_myfig = [];

            % expandfig request
            case {'expandfig'}

                % discard more than two size arguments
                if size(h) > 2
                    h = h(1:2);
                end

                % set new position
                fp = get(i_myfig.f, 'Position');
                fp(3:4) = fp(3:4) + h;
                set(i_myfig.f, 'Position', fp);

            % figsize request
            case {'figsize'}

                % get position
                fp = get(i_myfig.f, 'Position');

                % if size of argument is OK, set size
                if numel(h) == 2
                    fp(3:4) = h(:)';
                    set(i_myfig.f, 'Position', fp);
                end

                % otherwise get output
                if nargout > 0
                    varargout{1} = fp(3:4);
                end

            % figure request
            case {'figure'}
                varargout{1} = i_myfig.f;

            % image request
            case {'image'}
                varargout{1} = i_myfig.i;

            % newaxes request
            case {'newaxes'}

                % only valid for 1x4 array
                if numel(h) == 4

                    % add axes
                    i_myfig.h.a(end+1) = axes('Parent', i_myfig.f);

                    % set options
                    set(i_myfig.h.a(end),'Units',    'pixels', ...
                                         'Position', h(:)');

                    % return new object handle
                    varargout{1} = i_myfig.h.a(end);
                end

            % newfig request
            case {'newfig'}
                i_myfig = [];

            % text request
            case {'text'}
                varargout{1} = i_myfig.x;

            % visible request
            case {'visible'}

                % only valid for non empty h
                if isempty(h)
                    return;
                end

                % if first element true, otherwise ...
                if h(1)
                    set(i_myfig.f,'Visible','on');
                else
                    set(i_myfig.f,'Visible','off');
                end

            % non-request char argument
            otherwise

                % no internal structure as of yet
                if isempty(i_myfig)

                    % but handle given
                    if ~isempty(h)

                        % display on handle
                        dispslice(h,varargin{1:end});
                    end

                % with structure...
                else

                    % try first
                    try
                        get(i_myfig.i, 'Visible');

                    % and create new figure if needed
                    catch ne_eo;
                        neuroelf_lasterr(ne_eo);
                        dispslice('newfig');

                        % if now h was given
                        if ~isempty(h)

                            % display
                            dispslice(h, varargin{1:end});
                        end
                    end
                end

                % if still no struct
                if isempty(i_myfig)
                    return;
                end

                % otherwise set title (last resort for char argument)
                set(i_myfig.f, 'Name', img);

                % more arguments ?
                if nargin > 1

                    % don't accept empty arguments
                    if isempty(h)
                        return;
                    end

                    % if is char, only accept numbers
                    if ischar(h)
                        try
                            h = str2double(h);
                        catch ne_eo;
                            neuroelf_lasterr(ne_eo);
                            return;
                        end
                    end

                    % if one number is given, display zero square slice
                    if numel(h) == 1
                        dispslice(zeros(h, h));

                    % with two numbers, display zero rectangle slice
                    elseif numel(h) == 2
                        dispslice(zeros(h(:)'));

                    % otherwise this is the slice
                    else
                        dispslice(h);
                    end
                end
        end
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
    end

    % return now (with the first argument of type char!)
    return;
end

% reduce dimensions and force class
img = double(squeeze(img));

% only one argument, use default
if nargin < 2

    % no i_myfig yet
    if isempty(i_myfig)

        % create struct
        i_myfig.h = struct('a', [], 'i', [], 'l', [], 't', []);

        % make it the right size
        isz  = size(img);
        ifsz = max(isz);
        while ifsz > 512
            ifsz = ifsz / 2;
        end
        while ifsz < 256
            ifsz = ifsz * 2;
        end

        % create figure
        i_myfig.f = figure('Name',         'DispSlice', ...
                           'Tag',          'figDispSlice', ...
                           'MenuBar',      'none', ...
                           'NumberTitle',  'off', ...
                           'Position',     [64 64 15+ifsz 80+ifsz], ...
                           'DoubleBuffer', 'on');

        % add image
        i_myfig.i = image(zeros(isz));

        % set colormap to full grayscale
        colormap(repmat(0:1/255:1, [3, 1])');

        % set other optionos
        i_myfig.a = get(i_myfig.i, 'Parent');
        i_myfig.x = axes('Parent', i_myfig.f);
        i_myfig.t = text(0, 0, 'Image Info: ', 'Parent', i_myfig.x);

        % bring to front
        figure(i_myfig.f);

        % set current figure option
        set(i_myfig.f, 'HandleVisibility', 'off');

        % set object options
        set(i_myfig.a, 'Visible',  'off', ...
                       'Units',    'pixels', ...
                       'Position', [8 72 ifsz([1, 1])], ...
                       'HandleVisibility', 'off');
        set(i_myfig.t, 'Units',    'pixels', ...
                       'FontName', 'FixedWidth', ...
                       'Position', [8 40]);
        set(i_myfig.x, 'Visible',  'off', ...
                       'Units',    'pixels', ...
                       'Position', [6 0 ifsz(1) 64], ...
                       'HandleVisibility', 'off');

        % redraw (flush event queue)
        drawnow;
    end

    % try internal figure field
    try
        get(i_myfig.i, 'Visible');

    % recreate if necessary
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        i_myfig = [];
        dispslice('DispSlice', img);
        return;
    end

    % set h if needed
    h = i_myfig.i;
end

% reject empty or N-D data
if isempty(img) || ...
    length(size(img)) > 2
    return;
end

% nicely scale data to colormap if needed
mmm = minmaxmean(img, 1);
if nargin < 3 || ...
   ~isa(varargin{1}, 'double') || ...
    numel(varargin{1}) ~= 2
    imn = mmm(1);
    imx = mmm(2);

% or get min max from input
else
    imn = varargin{1}(1);
    imx = varargin{1}(2);
end

% give text as well ?
if h == i_myfig.i

    % get mean and std
    imm = mmm(3);
    ims = sqrt(mmm(6));

    % and set string
    set(i_myfig.t, 'String', ...
        sprintf(['Image Info:\n' ...
            'Rows:%4d, Cols:%4d, Pxls:%d\n' ...
            'Min: %9s, Max: %9s\n' ...
            'Mean:%9s, Std: %9s'], ...
        size(img, 1), size(img, 2), numel(img), ...
        sprintf('%.3g', imn), ...
        sprintf('%.3g', imx), ...
        sprintf('%.3g', imm), ...
        sprintf('%.3g', ims)));
end
img = uint8(floor(255.99 * (img - imn) / (realmin + (imx - imn))));
set(h, 'CData', img);
drawnow;
