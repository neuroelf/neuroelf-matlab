function xo = tom_ExtractSpots(xo, coordfile, opts)
% TOM::ExtractSpots  - extract spots from a file
%
% FORMAT:       tom.ExtractSpots(coordfile [, opts]);
%
% Input fields:
%
%       coordfile   JSON filename with coordinates
%       opts        1x1 struct with optional settings
%       .cutsize    optional cut size (default: 512 pixels)
%       .filter     filter expression, default: '$status==0'
%       .jpgqual    JPG quality (default: 90)
%       .marksize   optional mark size (default: 0 = no marking)
%       .outdir     output folder (default: 'spots' in TOM folder)
%
% No output fields.
%
% Using ddeblank.

% Version:  v1.1
% Build:    21102012
% Date:     Oct-20 2021, 12:50 PM EST
% Author:   Jochen Weber, NeuroElf.net, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2021, Jochen Weber
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

% check arguments
if nargin < 1 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'tom')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
t = xo.C;
tf = xo.F;
td = fileparts(tf);
if nargin < 2 || ~ischar(coordfile) || isempty(coordfile) || exist(coordfile(:)', 'file') ~= 2
    af = ne_methods.findfiles(td, 'analysis', '-D');
    cf = {};
    if numel(af) == 1
        cf = ne_methods.findfiles(af{1}, 'lesion*.json');
        if ~isempty(cf)
            coordfile = cf{end};
        end
    end
    if isempty(cf)
        error('neuroelf:xff:badArgument', 'Missing lesion analysis file.');
    end
end

% try loading coordfile
try
    coords = jsondecode(ne_methods.asciiread(coordfile(:)'));
    if ~isfield(coords, 'root') || ~isstruct(coords.root) || ...
        numel(coords.root) ~= 1 || ~isfield(coords.root, 'children')
        error('neuroelf:xff:badArgument', 'Invalid lesion analysis file.');
    end
catch ne_eo;
    rethrow(ne_eo);
end
cc = coords.root.children;
cf = fieldnames(cc);

% check options
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'cutsize') || ~isa(opts.cutsize, 'double') || numel(opts.cutsize) ~= 1 || ...
    isinf(opts.cutsize) || isnan(opts.cutsize) || opts.cutsize < 0
    opts.cutsize = 512;
else
    opts.cutsize = min(2048, max(64, round(opts.cutsize)));
end
if ~isfield(opts, 'filter') || ~ischar(opts.filter)
    opts.filter = '$status==0';
else
    opts.filter = opts.filter(:)';
end
if ~isfield(opts, 'jpgqual') || ~isa(opts.jpgqual, 'double') || numel(opts.jpgqual) ~= 1 || ...
    isinf(opts.jpgqual) || isnan(opts.jpgqual)
    opts.jpgqual = 90;
else
    opts.jpgqual = min(100, max(25, round(opts.jpgqual)));
end
if ~isfield(opts, 'marksize') || ~isa(opts.marksize, 'double') || numel(opts.marksize) ~= 1 || ...
    isinf(opts.marksize) || isnan(opts.marksize) || opts.marksize < 0
    opts.marksize = 0;
elseif opts.marksize > 0
    opts.marksize = min(ceil(0.5 * opts.cutsize), max(32, round(opts.marksize)));
end
if ~isfield(opts, 'outdir') || ~ischar(opts.outdir) || isempty(opts.outdir)
    opts.outdir = [td filesep 'spots'];
else
    opts.outdir = opts.outdir(:)';
    if opts.outdir(1) ~= filesep
        opts.outdir = ['./' opts.outdir];
    end
end
sf = opts.outdir;
if exist(sf, 'dir') ~= 7
    try
        ne_methods.mkadir(sf);
        if exist(sf, 'dir') ~= 7
            error('neuroelf:xff:mkadirError', 'Could not create spots folder.');
        end
    catch ne_eo;
        rethrow(ne_eo);
    end
end

% parse filter expression
cond = ne_methods.ddeblank(opts.filter);

% no condition, use full selection
if isempty(cond)
    sel = true(numel(cc), 1);
    
% otherwise
else

    % format cond as expression
    cond = ['(' cond ')'];

    % replace column names (==)
    cregx = regexp(cond, '(\$[a-zA-Z][a-zA-Z_0-9]*\s*==\s*''[^'']+'')', 'tokens');
    while ~isempty(cregx) && ~isempty(cregx{1})
        cregp = regexp(cregx{1}{1}, '^\$([a-zA-Z][a-zA-Z_0-9]*)\s*==\s*''([^'']+)''$', 'tokens');
        if ~iscell(cregp) || numel(cregp) ~= 1 || ~iscell(cregp{1}) || numel(cregp{1}) ~= 2
            error('neuroelf:xff:badArgument', 'Invalid conditional statement.');
        end
        if any(cregp{1}{2} == '*')
            cond = strrep(cond, cregx{1}{1}, sprintf( ...
                '~cellfun(''isempty'', regexpi({cc.%s}, ''%s''))', cregp{1}{:}));
        else
            cond = strrep(cond, cregx{1}{1}, sprintf( ...
                'strcmpi({cc.%s}, ''%s'')', cregp{1}{:}));
        end
        cregx = regexp(cond, '(\$[a-zA-Z][a-zA-Z_0-9]*\s*==\s*''[^'']+'')', 'tokens');
    end
    cregx = regexp(cond, '(\$[a-zA-Z][a-zA-Z_0-9]*\s*~=\s*''[^'']+'')', 'tokens');
    while ~isempty(cregx) && ~isempty(cregx{1})
        cregp = regexp(cregx{1}{1}, '^\$([a-zA-Z][a-zA-Z_0-9]*)\s*~=\s*''([^'']+)''$', 'tokens');
        if ~iscell(cregp) || numel(cregp) ~= 1 || ~iscell(cregp{1}) || numel(cregp{1}) ~= 2
            error('neuroelf:xff:badArgument', 'Invalid conditional statement.');
        end
        if any(cregp{1}{2} == '*')
            cond = strrep(cond, cregx{1}{1}, sprintf( ...
                'cellfun(''isempty'', regexpi({cc.%s}, ''%s''))', cregp{1}{:}));
        else
            cond = strrep(cond, cregx{1}{1}, sprintf( ...
                '~strcmpi({cc.%s}, ''%s'')', cregp{1}{:}));
        end
        cregx = regexp(cond, '(\$[a-zA-Z][a-zA-Z_0-9]*\s*~=\s*''[^'']+'')', 'tokens');
    end
    for op = {'==', '>=', '<=', '~=', '>', '<'}
        cregx = regexp(cond, ['(\$[a-zA-Z][a-zA-Z_0-9]*\s*' op{1} '\s*\-?[0-9\.\+\-eE]+)'], 'tokens');
        while ~isempty(cregx) && ~isempty(cregx{1})
            cregp = regexp(cregx{1}{1}, ['^\$([a-zA-Z][a-zA-Z_0-9]*)\s*' op{1} '\s*(\-?[0-9\.\+\-eE]+)$'], 'tokens');
            if ~iscell(cregp) || numel(cregp) ~= 1 || ~iscell(cregp{1}) || numel(cregp{1}) ~= 2
                error('neuroelf:xff:badArgument', 'Invalid conditional statement.');
            end
            cond = strrep(cond, cregx{1}{1}, sprintf(['[cc.%s] ' op{1} ' %s'], cregp{1}{:}));
            cregx = regexp(cond, ['(\$[a-zA-Z][a-zA-Z_0-9]*\s*' op{1} '\s*\-?[0-9\.\+\-eE]+)'], 'tokens');
        end
    end

    % then parse condition
    try
        sel = eval(cond);
    catch xfferror
        error('neuroelf:xff:badArgument', 'Bad condition given: %s.', xfferror.message);
    end
end
cc = cc(sel);

% iterate over selected spots
for c = 1:numel(cc)
    try
        msel = tom_MarkSpot(xo, [cc(c).x, cc(c).y, cc(c).z]);
    catch xfferror
        warning('neuroelf:xff:lookupError', 'Error locating spot at (%g/%g/%g): %s.', ...
            cc(c).x, cc(c).y, cc(c).z, xfferror.message);
        continue;
    end
    cut = tom_ExtractSpot(xo, msel, opts.cutsize, opts.marksize);
    imwrite(cut, sprintf('%s/%s.jpg', sf, cc(c).id), 'Quality', opts.jpgqual);
end
