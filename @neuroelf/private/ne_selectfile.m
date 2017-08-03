% FUNCTION ne_selectfile: populate Select current object menu
function ne_selectfile(varargin)

% Version:  v1.1
% Build:    16050816
% Date:     May-08 2016, 4:49 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2016, Jochen Weber
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

% global variable
global ne_gcfg;

% get handles for menu, voi, and voi index object
ch = ne_gcfg.h.Menu.SelectFile;
mh = ch.MLHandle;
scvarh = ne_gcfg.h.Scenery;
slvarh = ne_gcfg.h.SliceVar;
scvar = scvarh.UserData;
scidx = scvarh.Value;
slvar = slvarh.UserData;
slidx = slvarh.Value;

% for no further arguments
if nargin < 3 || ~ischar(varargin{3}) || isempty(varargin{3})

    % remove all children first
    delete(ch.Children);

    % only once at a time
    if ne_gcfg.c.incb
        return;
    end
    ne_gcfg.c.incb = true;

    % no or empty stvar
    if isempty(slvar) && isempty(scvar)

        % simply inform user...
        uimenu(mh, 'Label', 'No objects loaded', 'Enable', 'off');
        ne_gcfg.c.incb = false;
        return;
    end

    % fill menu
    if ~isempty(slvar)
        slvarr = false(size(slvar, 1), 1);
        uimenu(mh, 'Label', 'Slicing objects', 'Enable', 'off');
        for vc = 1:size(slvar, 1)
            slv = slvar{vc, 4};
            if numel(slv) ~= 1 || ~isxff(slv, true)
                slvarr(vc) = true;
                continue;
            end
            slvf = slv.FilenameOnDisk(true);
            if isempty(slvf)
                slvf = sprintf('<%s, xff #%d>', upper(slv.Filetype), slv.Filenumber);
            else
                [slvp, slvf] = fileparts(slvf);
            end
            mmh = uimenu(mh, 'Label', ['    ' slvf], 'Callback', {@ne_setcvar, vc});
            if vc == slidx
                set(mmh, 'Checked', 'on');
            end
        end
        if any(slvarr)
            slvart = slvarh.String;
            if ~iscell(slvart)
                slvart = cellstr(slvart);
            end
            slvart(slvarr) = [];
            slvar(slvarr, :) = [];
            slvarh.UserData = slvar;
            if isempty(slvart)
                slvart = {'empty'};
            end
            slvarh.String = slvart;
            if slvarr(slidx)
                slvarh.Value = 1;
            end
            try
                ne_gcfg.c.incb = false;
                ne_setcvar;
            catch ne_eo;
                ne_gcfg.c.lasterr = ne_eo;
            end
            ne_gcfg.c.incb = true;
        end
    else
        uimenu(mh, 'Label', 'No slicing objects loaded', 'Enable', 'off');
    end

    % first add tools entries
    if ~isempty(scvar)
        uimenu(mh, 'Label', 'Scenery objects', 'Enable', 'off', 'Separator', 'on');
        for vc = 1:size(scvar, 1)
            scv = scvar{vc, 4};
            if numel(scv) ~= 1 || ~isxff(scv)
                continue;
            end
            scvf = scv.FilenameOnDisk(true);
            if isempty(scvf)
                scvf = sprintf('<SRF #%d>', scv.Filenumber);
            else
                [scvp, scvf] = fileparts(scvf);
            end
            mmh = uimenu(mh, 'Label', ['    ' scvf], ...
                'Callback', {@ne_selectfile, 'scenery', vc});
            if any(scidx == vc)
                set(mmh, 'Checked', 'on');
            end
        end
    else
        uimenu(mh, 'Label', 'No surface objects loaded', 'Enable', 'off', 'Separator', 'on');
    end

    % remove only-once flag
    ne_gcfg.c.incb = false;
    return;
end

% only once at a time
if ne_gcfg.c.incb
    return;
end
ne_gcfg.c.incb = true;

% action
action = lower(varargin{3}(:)');
switch (action)

    % select a different cluster
    case {'scenery'}

        % valid selection
        if nargin < 4 || ~isa(varargin{4}, 'double') || numel(varargin{4}) ~= 1 || ...
            isinf(varargin{4}) || isnan(varargin{4}) || varargin{4} < 1 || ...
            varargin{4} > size(scvar, 1) || varargin{4} ~= fix(varargin{4})
            ne_gcfg.c.incb = false;
            return;
        end

        % flip flag
        if numel(scidx) > 1 && any(scidx == varargin{4})
            scidx(scidx == varargin{4}) = [];
        else
            scidx = sort([scidx(:); varargin{4}]);
        end

        % set new scenery content
        scvarh.Value = scidx;
        ne_gcfg.c.incb = false;
        ne_setsurfpos(0, 0, 1);
        ne_showpage(0, 0, 3);
        ne_setsurfpos;
        drawnow;
end

% remove only-once flag
ne_gcfg.c.incb = false;
