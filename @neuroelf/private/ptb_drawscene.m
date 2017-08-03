function sc = ptb_drawscene(ptbvars, scene)
% ptb_drawscene  - draw scene into buffer
%
% FORMAT:       [sc] = ptb_drawscene(ptbvars, scene)
%
% Input fields:
%
%       ptbvars     1x1 struct created by ptb_initscreen
%       scene       1xC string representing the scene
%
% Output fields:
%
%       sc          if requested, Cx1 cell array with 1xF cell arrays
%                   to use in feval calls to re-draw the scene
%
% Note: does NOT include the Screen('Flip') call!
%
% Note: the scene string written as follows:
%
% <TAB-CHAR><ELEMENT-CHAR>
% Type<TAB>Content<TAB>X<TAB>Y<TAB>[X2<TAB>Y2<TAB>]COLOR<TAB><ELEMENT>
% ...

% Version:  v0.9d
% Build:    14061709
% Date:     Jun-17 2014, 9:50 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2012, 2014, Jochen Weber
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

% allow special case
if nargin == 1 && ...
    iscell(ptbvars)

    % attempt to draw
    for c = 1:numel(ptbvars)
        if numel(ptbvars{c}) > 1 && ...
            iscell(ptbvars{c}) && ...
            isa(ptbvars{c}{1}, 'function_handle')
            try
                feval(ptbvars{c}{:});
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
            end
        end
    end

    % return early
    return;
end

% argument check
if nargin < 2 || ...
   ~isstruct(ptbvars) || ...
    numel(ptbvars) ~= 1 || ...
   ~isfield(ptbvars, 'screen') || ...
   ~isfield(ptbvars, 'screenrect') || ...
   ~isfield(ptbvars, 'texturenames') || ...
   ~isfield(ptbvars, 'textures') || ...
   ~isfield(ptbvars, 'texturesizes') || ...
   ~isa(ptbvars.screen, 'double') || ...
    numel(ptbvars.screen) ~= 1 || ...
   ~isa(ptbvars.screenrect, 'double') || ...
   ~isequal(size(ptbvars.screenrect), [1, 4]) || ...
   ~isstruct(ptbvars.textures) || ...
    numel(ptbvars.textures) ~= 1 || ...
   ~isstruct(ptbvars.texturesizes) || ...
    numel(ptbvars.texturesizes) ~= 1 || ...
   ~iscell(ptbvars.texturenames) || ...
    size(ptbvars.texturenames, 2) ~= 2 || ...
   ~ischar(scene) || ...
    numel(scene) < 2
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid of missing argument.' ...
    );
end

% get inputs
s = ptbvars.screen;
r = ptbvars.screenrect;
t = ptbvars.textures;
tn = ptbvars.texturenames;
ts = ptbvars.texturesizes;

% function handle
sfh = @Screen;

% get type and element characters
tchar = scene(1);
echar = scene(2);
scene = splittocellc(scene(3:end), echar);

% remove empty elements
scene(cellfun('isempty', scene)) = [];

% empty scene
if isempty(scene)

    % draw black dot in [0, 0, 0, 0]
    Screen('FillRect', s, [0, 0, 0, 0], [0, 0, 1, 1]);

    % output?
    if nargout > 0
        sc = {{sfh, 'FillRect', s, [0, 0, 0, 0], [0, 0, 1, 1]}};
    end

    % return
    return;
end

% valid elements
elements = { ...
    'DrawText'; ...
    'DrawTexture'; ...
    'FillRect'; ...
    'FrameRect'};

% create output
sc = scene(:);

% iterate over scene elements
sk = false(numel(sc), 1);
for c = 1:numel(sk)

    % split scene
    element = splittocellc(scene{c}, tchar);

    % only valid elements!
    if numel(elements) < 2 || ...
       ~any(strcmpi(element{1}, elements))
        continue;
    end

    % assume things are off
    fargs = {};

    % depending on element
    switch (lower(element{1}))

        % draw text
        case {'drawtext'}

            % we need the text and the X/Y positions
            if numel(element) < 4
                continue;
            end
            text = strrep(element{2}, '<crlf>', char(10));
            x = element{3};
            y = element{4};

            % test
            if isempty(text) || ...
                isempty(regexpi(x, '^(center|\d+)$')) || ...
                isempty(regexpi(y, '^(center|\d+)$'))
                continue;
            end

            % convert
            if x(1) <= 57
                x = str2double(x);
            end
            if y(1) <= 57
                y = str2double(y);
            end

            % color?
            if numel(element) > 4 && ...
               ~isempty(regexpi(element{5}, '^\s*\d+\s*,\s*\d+\s*,\s*\d+\s*$'))
                color = u8str2double(element{5}, 1, 3);
            else
                color = [255, 255, 255];
            end
            color = round(limitrangec(color, 0, 255, 0));

            % function arguments
            fargs = {@DrawFormattedText, s, text, x, y, color};

        % draw texture
        case {'drawtexture'}

            % at least texture and X/Y
            if numel(element) < 4
                continue;
            end
            texture = element{2};
            x = element{3};
            y = element{4};

            % test
            if isempty(texture) || ...
                isempty(regexpi(x, '^(center|\d+)$')) || ...
                isempty(regexpi(y, '^(center|\d+)$'))
                continue;
            end

            % texture must be loaded!
            if isfield(t, texture)
                tid = t.(texture);
                tsz = ts.(texture);
            else
                tlup = findfirst(strcmpi(tn(:, 1), texture));
                if isempty(tlup)
                    continue;
                end
                tid = tn{tlup, 2};
                tsz = ts.(tid);
                tid = t.(tid);
            end

            % convert
            if x(1) <= 57
                x = str2double(x);
            else
                x = round(0.5 * (r(1) + r(3) - tsz(1)));
            end
            if y(1) <= 57
                y = str2double(y);
            else
                y = round(0.5 * (r(2) + r(4) - tsz(2)));
            end

            % function arguments
            fargs = {sfh, 'DrawTexture', s, tid, [], [x, y, x + tsz(1), y + tsz(2)]};

        % fill rectangle
        case {'fillrect'}

            % exactly color and rect
            if numel(element) ~= 3
                continue;
            end
            color = element{2};
            rect = element{3};

            % test
            if isempty(regexpi(color, '^\s*\d+(,|\s+)\s*\d+(,|\s+)\s*\d+\s*$')) || ...
                isempty(regexpi(rect, '^\s*\d+(,|\s+)\s*\d+(,|\s+)\s*\d+(,|\s+)\s*\d+\s*$'))
                continue;
            end
            color = limitrangec(u8str2double(color, 1, 3), 0, 255, 0);
            rect = u8str2double(rect, 1, 4);
            if any(rect(3:4) < rect(1:2))
                continue;
            end

            % function arguments
            fargs = {sfh, 'FillRect', s, color, rect};

        % fill rectangle
        case {'framerect'}

            % exactly color and rect
            if numel(element) ~= 4
                continue;
            end
            color = element{2};
            rect = element{3};
            penwidth = element{4};

            % test
            if isempty(regexpi(color, '^\s*\d+(,|\s+)\s*\d+(,|\s+)\s*\d+\s*$')) || ...
                isempty(regexpi(rect, '^\s*\d+(,|\s+)\s*\d+(,|\s+)\s*\d+(,|\s+)\s*\d+\s*$')) || ...
                isempty(regexpi(penwidth, '^\s*[1-9][0-9]*\s*$'))
                continue;
            end
            color = limitrangec(u8str2double(color, 1, 3), 0, 255, 0);
            rect = u8str2double(rect, 1, 4);
            penwidth = u8str2double(penwidth, 1, 1);
            if any(rect(3:4) < rect(1:2)) || ...
                penwidth > 16
                continue;
            end

            % function arguments
            fargs = {sfh, 'FrameRect', s, color, rect, penwidth};
    end

    % do something
    if ~isempty(fargs)
        try
            feval(fargs{:});
            sk(c) = true;
            if nargout > 0
                sc{c} = fargs;
            end
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
        end
    end
end

% keep the ones we used
sc = sc(sk);
