function tidx = pfs_LocateProperty(xo, pname, setdir)
% PFS::LocateProperty  - locate a property within a POIFS object
%
% FORMAT:       pidx = pfs.LocateProperty(pname)
%
% Input fields:
%
%       pname       property name
%
% Output fields:
%
%       pidx        property index (0 if not found)
%
% Using: findfirst, splittocellc.

% Version:  v1.1
% Build:    16021018
% Date:     Feb-10 2016, 6:10 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
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

% neuroelf library
global ne_methods;
findfirst = ne_methods.findfirst;

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'pfs') || ...
   ~ischar(pname) || isempty(pname) || (~any(pname(:) == '/') && numel(pname) > 31) || ...
    numel(pname) > 511
    error('neuroelf:xff:badArgument', 'Invalid argument in call.');
end
if nargin < 3 || ~islogical(setdir) || numel(setdir) ~= 1
    setdir = false;
end
bc = xo.C;
if isempty(bc.Properties)
    pfs_GetAllProperties(xo);
    bc = xo.C;
end
pname = pname(:)';

% check that CWD is good
if numel(bc.CWD) ~= 1 || ~isa(bc.CWD, 'double') || isnan(bc.CWD) || bc.CWD < 1 || ...
    bc.CWD > numel(bc.Properties) || ~any(bc.Properties(round(bc.CWD)).Type == [1, 5])
    bc.CWD = 1;
    xo.C = bc;
    error('neuroelf:xff:internalError', 'Invalid CWD in POIFS object set.');
end
ocwd = bc.CWD;

% sub directory from current
if any(pname == '/')
    pchain = ne_methods.splittocellc(pname, '/');
    if isempty(pchain{1})
        pchain(1) = [];
        if bc.CWD ~= 1
            bc.CWD = 1;
            xo.C = bc;
        end
    end
    while numel(pchain) > 1
        try
            pfs_LocateProperty(xo, pchain{1}, true);
            pchain(1) = [];
        catch xfferror
            bc.CWD = ocwd;
            xo.C = bc;
            rethrow(xfferror);
        end
    end
    if ~isempty(pchain)
        pname = pchain{1};
        bc = xo.C;
    else
        tidx = 1;
        return;
    end
end

% now try to locate the property
pnlen = numel(pname);
ps = bc.Properties;
tidx = ps(bc.CWD).FirstChild;
if numel(tidx) ~= 1 || ~isa(tidx, 'double') || isnan(tidx) || tidx < 1 || tidx > numel(bc.Properties)
    error('neuroelf:xff:internalError', 'Invalid FirstChild property for directory.');
end

% search for child
while ~strcmp(pname, ps(tidx).Name)
    if numel(ps(tidx).Name) < pnlen
        tidx = ps(tidx).Next;
    elseif numel(ps(tidx).Name) > pnlen
        tidx = ps(tidx).Previous;
    else
        pdiff = double(ps(tidx).Name) - double(pname);
        pdfs = pdiff(findfirst(pdiff ~= 0));
        if pdfs < 0
            tidx = ps(tidx).Next;
        else
            tidx = ps(tidx).Previous;
        end
    end
    if numel(tidx) ~= 1 || ~isa(tidx, 'double') || isnan(tidx) || tidx > numel(ps)
        error('neuroelf:xff:internalError', 'Invalid Previous/Next property detected.');
    end
    if tidx < 1
        if ~setdir
            tidx = 0;
            return;
        else
            error('neuroelf:xff:dirNotFound', 'Cannot set CWD; property not found.');
        end
    end
end
if setdir
    bc.CWD = tidx;
    xo.C = bc;
end
