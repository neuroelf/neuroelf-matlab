% FUNCTION cv_removefromlist: remove one variable from the workspace
function cv_removefromlist(h, v, emptystring)

% Version:  v1.1
% Build:    16012418
% Date:     Jan-24 2016, 6:44 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2014, 2016, Jochen Weber
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
if nargin < 3
    emptystring = 'empty';
end

% get current list
hc = get(h.MLHandle);
s = hc.String;
if ~iscell(s)
    s = cellstr(s);
end
s = s(:);

% default is to keep control on
e = 'on';

% if control is already empty
if numel(s) == 1 && ...
    strcmp(s{1}, emptystring)

    % then do nothing
    return;
end

% otherwise we need current data
u = hc.UserData;
x = hc.Value;

% iterate over the userdata
for c = size(u, 1):-1:1

    % equally create a struct
    ts = u{c, 4};

    % and compare xff handles, upon match
    if v == ts

        % remove from list
        u(c, :) = [];
        s(c) = [];
        if any(x == c)
            xc = find(x == c);
            x(xc) = [];
            if xc <= numel(x)
                x(xc:end) = x(xc:end) - 1;
            end
        else
            xc = (x > c);
            x(xc) = x(xc) - 1;
        end

        % but if nothing remains
        if isempty(s)

            % set to "empty" again
            s = {emptystring};
            if strcmpi(hc.Style, 'popupmenu') || ...
                (h.Max - h.Min) < 2
                x = 1;
            else
                x = [];
            end
            e = 'off';
        elseif isempty(x) && ...
            strcmpi(hc.Style, 'popupmenu') || ...
            (h.Max - h.Min) < 2
            x = max(1, c - 1);
        end

        % then set to control
        if strcmpi(hc.Style, 'listbox')
            set(h.MLHandle, ...
                'String',     s, ...
                'Value',      x, ...
                'ListboxTop', max(1, min(hc.ListboxTop, numel(s))), ...
                'UserData',   u, ...
                'Enable',     e);
        else
            set(h.MLHandle, ...
                'String',   s, ...
                'Value',    x, ...
                'UserData', u, ...
                'Enable',   e);
        end

        % and return early
        v.SetHandle('ShownInGUI', false);
        return;
    end
end
