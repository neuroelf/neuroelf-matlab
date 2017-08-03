function [f, p, l, t] = tar_FilesByContent(xo, n, opts)
% TAR::FilesByContent  - look for files by content
%
% FORMAT:       [files, pos, lengths, tio] = tar.FilesByContent(c [, opts])
%
% Input fields:
%
%       c           byte sequence (1xN uint8/char) to select files
%       opts        1x1 struct with optional settings
%        .icase     1x1 logical flag, ignore case (default: true)
%        .progress  either {false} or a 1x1 xfigure::progress or xprogress
%        .transio   1x1 logical flag, return cell of transios (false)
%        .within    1x2 cell of tokens, only look within boundaries
%                   both of these MUST return only a single occurrence
%
% Output fields:
%
%       files       filenames that match
%       pos         position within the TAR file
%       lengths     file lengths
%       tio         cell array with transio objects

% Version:  v1.1
% Build:    16060323
% Date:     Jun-03 2016, 11:13 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2016, Jochen Weber
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

% only valid for single file
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'tar') || ...
   (~isa(n, 'uint8') && ~ischar(n) && ~iscell(n)) || isempty(n)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
if ~iscell(n)
    n = {n};
end
n = n(:)';

% options
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'icase') || ~islogical(opts.icase) || numel(opts.icase) ~= 1
    opts.icase = ischar(n{1});
elseif ~ischar(n{1})
    opts.icase = false;
end
icase = opts.icase;
if ~isfield(opts, 'progress') || numel(opts.progress) ~= 1 || ...
   (~islogical(opts.progress) && ~isa(opts.progress, 'xprogress') && ~isxfigure(opts.progress, true))
    opts.progress = true;
end
if ~isfield(opts, 'transio') || ~islogical(opts.transio) || numel(opts.transio) ~= 1
    opts.transio = false;
end
if ~isfield(opts, 'within') || ~iscell(opts.within) || numel(opts.within) ~= 2 || ...
    isempty(opts.within{1}) || (~ischar(opts.within{1}) && ~isa(opts.within{1}, 'uint8')) || ...
    isempty(opts.within{2}) || (~ischar(opts.within{2}) && ~isa(opts.within{2}, 'uint8'))
    win = {};
else
    win = opts.within(:)';
    if ischar(win{1})
        if icase
            win1 = uint8(upper(win{1}(:)'));
        else
            win1 = uint8(win{1}(:)');
        end
    else
        win1 = win{1}(:)';
    end
    if ischar(win{2})
        if icase
            win2 = uint8(upper(win{2}(:)'));
        else
            win2 = uint8(win{2}(:)');
        end
    else
        win2 = win{2}(:)';
    end
    w1n1 = win1(1);
    w1nl = numel(win1);
    w1onl = ones(1, w1nl);
    w1cnl = 0:(w1nl-1);
    w2n1 = win2(1);
    w2nl = numel(win2);
    w2onl = ones(1, w2nl);
    w2cnl = 0:(w2nl-1);
end
n1 = zeros(numel(n), 1);
nl = zeros(numel(n), 1);
onl = cell(numel(n), 1);
cnl = cell(numel(n), 1);
for nc = 1:numel(n)
    if isempty(n{nc}) || (~ischar(n{nc}) && ~isa(n{nc}, 'uint8'))
        error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
    end
    n{nc} = n{nc}(:)';
    if ischar(n{nc})
        if icase
            n{nc} = upper(n{nc});
        end
        n{nc} = uint8(n{nc});
    end
    n1(nc) = n{nc}(1);
    nl(nc) = numel(n{nc});
    onl{nc} = ones(1, nl(nc));
    cnl{nc} = 0:(nl(nc)-1);
end
tnc = nc;

% determine progress bar capabilities
try
    closepbar = false;
    pbar = [];
    nextprog = now + 1 / 86400;
    if islogical(opts.progress) && opts.progress
        pbar = xprogress;
        xprogress(pbar, 'setposition', [80, 200, 640, 36]);
        xprogress(pbar, 'settitle', 'Content-searching TAR archive...');
        xprogress(pbar, 0, 'Searching in files...', 'visible', 0, 1);
        closepbar = true;
    elseif ~islogical(opts.progress)
        pbar = opts.progress;
        pbarvis = pbar.Visible;
        pbar.Progress(0, 'Searching in files...');
        pbar.Visible = 'on';
    end
catch xfferror
    pbar = [];
    neuroelf_lasterr(xfferror);
end

% match
bc = xo.C;
fid = fopen(xo.F, 'r');
if fid < 1
    error('neuroelf:xff:fileOpenError', 'Error opening TAR file for reading.');
end
m = false(numel(bc.Name), 1);
p = bc.ContPos - 1;
l = bc.ContLen;
t = bc.Type;
for fc = 1:numel(m)
    if t(fc) ~= '0'
        continue;
    end
    fseek(fid, p(fc), -1);
    lc = l(fc);
    if icase
        fcont = uint8(upper(fread(fid, [1, lc], 'uint8=>char')));
    else
        fcont = fread(fid, [1, lc], 'uint8=>uint8');
    end
    if ~isempty(pbar) && now > nextprog
        [nullp, fname, fext] = fileparts(bc.Name{fc});
        pbar.Progress(fc / numel(m), sprintf('Searching in %s%s...', fname, fext));
        nextprog = now + 1 / 86400;
    end
    if ~isempty(win)
        w1fm1 = find(fcont == w1n1);
        if isempty(w1fm1)
            continue;
        end
        w2fm1 = find(fcont == w2n1);
        if isempty(w2fm1)
            continue;
        end
        w1fme = w1fm1 + w1nl - 1;
        w1fetl = (w1fme > lc);
        if any(w1fetl)
            w1fm1(w1fetl) = [];
            if isempty(w1fm1)
                continue;
            end
        end
        w2fme = w2fm1 + w2nl - 1;
        w2fetl = (w2fme > lc);
        if any(w2fetl)
            w2fm1(w2fetl) = [];
            if isempty(w2fm1)
                continue;
            end
        end
        w1nm = numel(w1fm1);
        w2nm = numel(w2fm1);
        w1fp = find(all(fcont(w1fm1(:) * w1onl + ones(w1nm, 1) * w1cnl) == repmat(win1, w1nm, 1), 2));
        w2fp = find(all(fcont(w2fm1(:) * w2onl + ones(w2nm, 1) * w2cnl) == repmat(win2, w2nm, 1), 2));
        if numel(w1fp) ~= 1 || numel(w2fp) ~= 1 || w1fp >= w2fp
            continue;
        end
        fcont = fcont(w1fm1(w1fp)+w1nl:w2fm1(w2fp)-1);
        lc = numel(fcont);
    end
    nmatch = true;
    for nc = 1:tnc
        fm1 = find(fcont == n1(nc));
        if isempty(fm1)
            nmatch = false;
            break;
        end
        fme = fm1 + nl(nc) - 1;
        fetl = (fme > lc);
        if any(fetl)
            fm1(fetl) = [];
            if isempty(fm1)
                nmatch = false;
                break;
            end
        end
        nm = numel(fm1);
        if ~any(all(fcont(fm1(:) * onl{nc} + ones(nm, 1) * cnl{nc}) == repmat(n{nc}, nm, 1), 2))
            nmatch = false;
            break;
        end
    end
    m(fc) = nmatch;
end
fclose(fid);
m = find(m);

% nothing found
if isempty(m)
    f = cell(0, 1);
    p = zeros(0, 1);
    l = zeros(0, 1);
    t = cell(0, 1);

% found
else

    % remove folders
    m(xo.C.Type(m) > '0') = [];

    % get files, positions and lengths
    f = xo.C.Name(m);
    p = xo.C.ContPos(m);
    l = xo.C.ContLen(m);
    if nargout > 3
        t = cell(numel(m), 1);

        % transios?
        if opts.transio
            tio_obj = struct(transio(xo.F, 'ieee-le', 'uint8', 0, [1, 512]));
            tio_pos = p - 1;
            for fc = 1:numel(m)
                tio_obj.IOOffset = tio_pos(fc);
                tio_obj.DataDims = [1, l(fc)];
                t{fc} = transio(0, 'makeobject', tio_obj);
            end
        end
    end
end

% reset progress
if ~isempty(pbar)
    if closepbar
        closebar(pbar);
    else
        pbar.Visible = pbarvis;
    end
end
