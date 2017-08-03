% FUNCTION ne_keypress_uic: handle keyboard of UICs
function ne_keypress_uic(src, ke, varargin)

% Version:  v0.9b
% Build:    11050712
% Date:     Apr-09 2011, 11:48 PM EST
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

% get Key and Modifier from keyboard event (see Matlab docu!)
kk = ke.Key;
mn = ke.Modifier;

% determine which modifiers are pressed
km = false(1, 4);
if ~isempty(mn)
    try
        km = [ ...
            any(strcmpi('alt', mn)), ...
            any(strcmpi('control', mn)), ...
            any(strcmpi('shift', mn)), ...
            any(strcmpi('command', mn))];
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
    end
end

% get type of UIC
uicprops = get(src);
if ~strcmpi(uicprops.Type, 'uicontrol')
    return;
end
uictype = lower(uicprops.Style);

% depending on type
switch (uictype)

    % edit fields
    case {'edit'}

        % no modifier
        if ~any(km)

            % what key
            switch (lower(kk))

                % down/up arrow key
                case {'downarrow', 'uparrow'}

                    % see if it's a number
                    if ischar(uicprops.String) && ...
                       ~isempty(uicprops.String) && ...
                       ~isempty(regexpi(uicprops.String, ...
                            '^\s*[+\-]?\d+(\.\d*)?(e[+\-]?\d+)?\s*$'))
                        uicnum = str2double(uicprops.String);

                        % de/increase number
                        if lower(kk(1)) == 'd'
                            set(src, 'String', sprintf('%g', uicnum - 1));
                        else
                            set(src, 'String', sprintf('%g', uicnum + 1));
                        end

                        % callback ?
                        if nargin > 3 && ...
                            isa(varargin{2}, 'function_handle')
                            feval(varargin{2:end});
                        end
                    end
            end
        end
end
