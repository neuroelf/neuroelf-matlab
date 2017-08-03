function ps = pfs_DirList(xo, pat, opts)
% PFS::DirList  - directory list
%
% FORMAT:       [d = ] pfs.DirList([pat [, opts]])
%
% Input fields:
%
%       pat         optional pattern
%       opts        optional settings
%        .show      show dir to console (default: true)
%        .subdirs   if true, also traverse into subdirs (default: false)
%
% Output fields:
%
%       d           if requested a Nx2 cell array with names and IDs
%
% Note: if no output is requested, the directory with be printed instead

% Version:  v1.1
% Build:    16021018
% Date:     Feb-10 2016, 6:15 PM EST
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

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'pfs')
    error('neuroelf:xff:badArgument', 'Invalid call to PFS method.');
end
if nargin < 2 || ~ischar(pat) || isempty(pat)
    pat = '';
else
    pat = pat(:)';
end
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'show') || numel(opts.show) ~= 1 || ~islogical(opts.show)
    dbs = dbstack('-completenames');
    while ~isempty(dbs)
        if ~isempty(strfind(lower(dbs(1).file), '/@xff'))
            dbs(1) = [];
        end
    end
    if isempty(dbs)
        opts.show = true;
        opts.return = false;
    else
        opts.show = false;
    end
end
if ~isfield(opts, 'subdirs') || numel(opts.subdirs) ~= 1 || ~islogical(opts.subdirs)
    opts.subdirs = false;
end
bc = xo.C;
if isempty(bc.Properties)
    pfs_GetAllProperties(xo);
    bc = xo.C;
end
cwd = bc.CWD;
if numel(cwd) ~= 1 || ~isa(cwd, 'double') || cwd < 1 || cwd > numel(bc.Properties)
    bc.CWD = 1;
    xo.C = bc;
    error('neuroelf:xff:internalError', 'Invalid CWD set.');
end

% traverse all next/previous we find
ps = bc.Properties;
li = zeros(1, 64);
li(1) = ps(cwd).FirstChild;
ni = li;
nf = 1;
nl = 1;
ln = ni(1:nf);
while ~isempty(ln)
    nf = 0;
    ni(:) = 0;
    for lc = 1:numel(ln)
        p = ps(ln(lc));
        if p.Previous > 0
            nf = nf + 1;
            ni(nf) = p.Previous;
        end
        if p.Next > 0
            nf = nf + 1;
            ni(nf) = p.Next;
        end
    end
    if (nl + nf) > numel(li)
        li(nl+nf+64) = 0;
    end
    ln = ni(1:nf);
    li(nl+1:nl+nf) = ln;
    nl = nl + nf;
end
li(nl+1:end) = [];

% create list of property names
ps = ps(li);

% add subdirectories as well
if opts.subdirs
    sd = find([ps(:).Type] == 1);
    if ~isempty(sd)
        subopts = opts;
        subopts.show = false;
        sdc = cell(numel(sd), 1);
        sdn = numel(li);
        for sc = 1:numel(sd)
            try
                bc.CWD = li(sd(sc));
                xo.C = bc;
                sdc{sc} = pfs_DirList(xo, pat, subopts);
                sdn = sdn + numel(sdc{sc});
                for ssc = 1:numel(sdc{sc})
                    sdc{sc}(ssc).Name = regexprep(sdc{sc}(ssc).Name, ...
                        '^(.*)$', sprintf('%s/$1', bc.Properties(li(sd(sc))).Name));
                end
            catch xfferror
                neuroelf_lasterr(xfferror);
            end
        end
        bc.CWD = cwd;
        xo.C = bc;
        pnn = ps(1, ones(1, sdn));
        ntc = 1;
        for sc = 1:numel(sd)
            if ~isempty(sdc{sc})
                pnn(ntc:ntc+numel(sdc{sc})-1) = sdc{sc};
                ntc = ntc + numel(sdc{sc});
            end
        end
        pnn(ntc:end) = ps;
        ps = pnn;
    end
end

% filter ?
if ~isempty(pat)
    pnkeep = false(numel(ps), 1);
    for nn = 1:numel(pnkeep)
        pnkeep(nn) = isempty(regexpi(ps(nn).Name, pat));
    end
    ps(pnkeep) = [];
end

% sort
[pnames, pns] = sort({ps(:).Name});
ps = ps(pns);

% show
if opts.show
    for nn = 1:numel(ps)
        p = ps(nn);
        if p.FirstChild > 0
            p.Name(end+1) = '/';
            sz = '  DIR  ';
        else
            switch (p.BlockAccess(1))
                case 'B'
                    if p.Size >= 1024000
                        szm = '%5.1fM';
                        szv = p.Size / 1048576;
                    elseif p.Size > 50000
                        szm = '%5.1fk';
                        szv = p.Size / 1024;
                    else
                        szm = '%5d ';
                        szv = p.Size;
                    end
                    sz = sprintf(['%6dB  ' szm], p.FirstBlock, szv);
                case 'S'
                    sz = sprintf('%6dS  %5d ', p.FirstBlock, p.Size);
                otherwise
                    sz = '??????-       0?';
            end
        end
        pnames{nn} = sprintf('  %-31s  %s', p.Name, sz);
    end
    ps = char(pnames(:));
end
