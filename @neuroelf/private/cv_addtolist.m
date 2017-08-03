% FUNCTION cv_addtolist: add a variable to the workspace (ne_gcfg.w)
function cv_addtolist(h, f, n, xt, emptystring)

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:16 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2014, 2016, Jochen Weber
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

% check argument number
if nargin < 4
    xt = {};
end
if nargin < 5
    emptystring = 'empty';
end

% check for empty filenames (unsaved objects)
[fn{1:3}] = fileparts(f.FilenameOnDisk(true));
if isempty(fn{2})

    % make them "untitled"
    fn{2} = sprintf('<untitled.%s>', f.Filetype);

    % and prepend extension with '.'
    fn{3} = ['.' f.Filetype];

% otherwise
else

    % get complete name
    fn{2} = [fn{2} fn{3}];
end

% get string and userdata from dropdown
s = h.String;
if ~iscell(s)
    s = cellstr(s);
end
s = s(:);
u = h.UserData;

% take care of first addition
if numel(s) == 1 && ...
    strcmp(s{1}, emptystring)
    s(1) = [];
end

% if already in list, do nothing
for uc = 1:size(u)
    if u{uc, 4} == f
        return;
    end
end

% add fields?
if ~isempty(xt)
    for xc = numel(xt):-1:1
        try
            if ischar(xt{xc})
                xfc = f.(xt{xc});
                if isnumeric(xfc)
                    if numel(xfc) == 1
                        if xfc == round(xfc)
                            xff = '%d';
                        else
                            xff = '%g';
                        end
                    elseif numel(xfc) < 4
                        xfc = sprintf('%g ', xfc);
                        xfc = sprintf('(%s)', xfc(1:end-1));
                        xt{xc} = sprintf('%s: %s', xt{xc}, xfc);
                        continue;
                    else
                        xt{xc} = sprintf('%s: (%d %s values)', xt{xc}, ...
                            numel(xfc), class(xfc));
                        continue;
                    end
                elseif ischar(xfc)
                    xff = '%s';
                    xfc = xfc(:)';
                else
                    xt(xc) = [];
                    continue;
                end
                xt{xc} = sprintf(['%s: ' xff], xt{xc}, xfc);
            elseif iscell(xt{xc}) && ...
                numel(xt{xc}) > 1 && ...
                ischar(xt{xc}{1})
                xt{xc} = sprintf(xt{xc}{:});
            end
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            xt(xc) = [];
        end
    end
    if ~isempty(xt)
        xt = sprintf('%s, ', xt{:});
        xt = sprintf('(%s)', xt(1:end-2));
    else
        xt = '';
    end
    fn{2} = [fn{2} ' ' xt];
end

% then add to list
s{end+1} = fn{2};
u(end + 1, :) = {s{end}, lower(fn{3}(2:end)), n, f};

% and set back to control
h.String = s;
h.UserData = u;

% then ensure that future calls work
f.SetHandle('ShownInGUI', true);
if strcmpi(h.Style, 'popupmenu')
    h.Value = numel(s);
elseif strcmpi(h.Style, 'listbox')
    h.ListboxTop = min(numel(s), h.ListboxTop);
end
h.Enable = 'on';
