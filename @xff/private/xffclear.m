function xffclear(l, cgui)
% xff::_clear  - clears objects' memory for given handles

% Version:  v1.1
% Build:    16051812
% Date:     May-18 2016, 12:29 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2015, 2016, Jochen Weber
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

% global storage
global xffsngl;
global ne_gcfg;

% nothing to be done
if isempty(xffsngl) || isempty(xffsngl.OBJS)
    return;
end

% otherwise try
try

    % for single ID, pack into cell
    if ischar(l)
        l = {l};
    end
    
    % get indices of objects to clear (and disallow ROOT)
    [rl{1:2}] = intersect(xffsngl.OBJS(:, 3), l);
    rl = rl{2};
    rl(rl < 2) = [];
    if isempty(rl)
        return;
    end
    
    % get objects
    clu = xffsngl.OBJS(rl, 4);

    % now process handles
    for c = 1:numel(rl)
        h = clu{c}.H;
        if isfield(h, 'CleanUp') && iscell(h.CleanUp) && ~isempty(h.CleanUp)
            clus = h.CleanUp;
            for cluc = 1:numel(clus)
                if ischar(clus{cluc}) && ~isempty(clus{cluc})
                    evalin('base', clus{cluc}(:)', '');
                elseif ishandle(clus{cluc}) || ...
                   (numel(clus{cluc}) == 1 && ~any(strcmpi(class(clus{cluc}), ...
                    {'int8', 'int16', 'int32', 'int64', 'uint8','uint16','uint32','uint64', ...
                     'char','single','double','struct', 'cell','function_handle'})))
                    try
                        delete(clus{cluc});
                    catch xfferror
                        neuroelf_lasterr(xfferror);
                    end
                end
            end
        end
        if isfield(h, 'ShownInGUI') && islogical(h.ShownInGUI) && ...
            numel(h.ShownInGUI) == 1 && h.ShownInGUI && (nargin < 2 || cgui) && ...
           ~isempty(ne_gcfg)
            neuroelf_gui('closefile', clu{c}, false);
        end
        if isfield(h, 'GZIPext') && ischar(h.GZIPext) && strcmpi(h.GZIPext, '.gz') && ...
            isfield(h, 'GZIPfile') && ischar(h.GZIPfile) && ~isempty(h.GZIPfile) && ...
            exist([h.GZIPfile h.GZIPext], 'file') == 2
            cluf = clu{c}.F;
            if ~isempty(cluf) && exist(cluf, 'file') == 2
                try
                    delete(cluf);
                catch xfferror
                    neuroelf_lasterr(xfferror);
                end
            end
        end
    end
    
    % and then remove from global storage
    xffsngl.OBJS(rl, :) = [];

catch xfferror
    error('neuroelf:xff:lookupError', ...
        'Error looking up/clearing objects: %s.', xfferror.message);
end
