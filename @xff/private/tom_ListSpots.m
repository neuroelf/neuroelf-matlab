function l = tom_ListSpots(xo, opts)
% TOM::ListSpots  - list spots from default JSON analysis file
%
% FORMAT:       tom.ExtractSpots([opts]);
%
% Input fields:
%
%       opts        1x1 struct with optional settings
%       .boxwidth   1x1 factor for box width (default: 1.25)
%       .delempty   1x1 boolean flag, remove empty/failed entries (true)
%       .filter     filter expression, default: '$status==0'
%       .outputs    list of fields to list, default:
%                   {'cx', 'cy', 'majorAxisMM', 'minorAxisMM', ...
%                    'deltaLBnorm', 'norm_border', 'norm_color', 'H', ...
%                    'nevi_confidence', 'location_simple'}
%
% Output fields:
%
%       l           CxF cell with list of selected spots
%
% Using findfiles, multimatch, splittocellc.

% Version:  v1.1
% Build:    22061513
% Date:     Jun-15 2022, 1:00 PM EST
% Author:   Jochen Weber, NeuroElf.net, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2022, Jochen Weber
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

% options
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'boxwidth') || ~isa(opts.boxwidth, 'double') || numel(opts.boxwidth) ~= 1 || ...
    isinf(opts.boxwidth) || isnan(opts.boxwidth) || opts.boxwidth < 0.7
    opts.boxwidth = 1.25;
elseif opts.boxwidth > 2.5
    opts.boxwidth = 2.5;
end
if ~isfield(opts, 'delempty') || ~islogical(opts.delempty) || numel(opts.delempty) ~= 1
    opts.delempty = true;
end
if ~isfield(opts, 'filter') || ~ischar(opts.filter)
    opts.filter = '$status==0';
else
    opts.filter = opts.filter(:)';
end
if ~isfield(opts, 'output') || ~iscell(opts.output)
    opts.output = {'cx'; 'cy'; 'majorAxisMM'; 'minorAxisMM'; 'deltaLBnorm'; ...
        'norm_border'; 'norm_color'; 'H'; 'nevi_confidence'; 'location_simple'};
else
    opts.output = opts.output(:);
    opts.output(~cellfun(@ischar, opts.output)) = [];
    opts.output = ne_methods.ddeblank(outs.output);
end
fo = {'json_index'; 'texture_id'; 'texture_width'; 'texture_height'; ...
    't_x1'; 't_x2'; 't_y1'; 't_y2'};
nf = numel(fo);
bc = xo.C;

% extract camera model names (for texture filenames)
f = bc.Field;
fn = {f.ContentName};
txs = find(strcmpi(fn, 'txtrjpg_') | strcmpi(fn, 'txtrjpga'));
md = bc.MetaData;
mn = {md.Name};
mi = find(strcmpi(mn(:), 'texcameramodelnames'));
if numel(mi) ~= 1
    warning('Cannot guarantee camera/image names.');
    l = [fo', opts.output'];
    return
else
    modelnames = ne_methods.splittocellc(md(mi).Content, char(0));
    if numel(modelnames) < numel(txs)
        for txc = (numel(modelnames)+1):numel(txs)
            modelnames{txc} = sprintf('texture%d', txc-1);
        end
    end
end
csz = struct;

% try to read JSON file
tom_folder = fileparts(xo.F);
a_folder = [tom_folder '/analysis'];
tx_folder = [tom_folder '/textures'];
use_textures = true;
if exist(tx_folder, 'dir') ~= 7
    warning('neureolf:missingData', 'Cannot read texture images.');
    use_textures = false;
    imsz = [1, 1, 3];
end
afiles = sort(ne_methods.findfiles(a_folder, 'lesion_data*.json', '-d1'));
if numel(afiles) < 1
    fprintf('No analysis folder/data found.');
    l = [fo', opts.output'];
    return
end
try
    tj = jsondecode(ne_methods.asciiread(afiles{end}));
    cc = tj.root.children;
catch ne_eo;
    rethrow(ne_eo);
end

% no condition, use full selection
cond = ne_methods.ddeblank(opts.filter);
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
sel = find(sel);
cc = cc(sel);
nc = numel(cc);
if ~isempty(opts.output)
    ofn = fieldnames(cc);
    ofm = ne_methods.multimatch(lower(opts.output), lower(ofn));
    opts.output = ofn(ofm(ofm > 0));
end
no = numel(opts.output);
l = cell(nc, nf + no);

% access object
crd = bc.VertexCoordinate;
tr = bc.TriangleVertex;
trac = bc.TexVertACoord;
tram = bc.TexVertAMap;
trm = bc.CornerTexVtxMap;

% iterate over children
for c = 1:nc
    
    % get specification
    icc = cc(c);
    ccrd = [icc.x, icc.y, icc.z];
    try
        [ms, vs] = tom_MarkSpot(xo, ccrd);
    catch ne_eo;
        warning('neuroelf:invalidData', 'Coordinate not close enough to any vertex.');
        l(c, :) = repmat({''}, 1, size(l, 2));
        continue;
    end
    
    % also get coordinates of triangle points
    trvtx = tr(vs(2), :);
    trcrd = crd(trvtx, :);
    txvtx = trm(trm(:, 1) == vs(2), 2:3);
    [~, txord] = sort(txvtx(:, 1));
    txvtx = txvtx(txord, 2);
    txcrd = double(trac(txvtx, :));
    txmap = tram(txvtx);
    if ~all(txmap == ms(1))
        warning('neuroelf:invalidData', 'Coordinate not close enough to any vertex.');
        l(c, :) = repmat({''}, 1, size(l, 2));
        continue;
    end
    trdst = sqrt(sum((trcrd([1,1], :) - trcrd(2:3, :)) .^ 2, 2));
    
    % compare index
    msm = modelnames{ms(1)};
    if ~isfield(csz, msm)
        if use_textures
            try
                im = imread([tx_folder '/' msm '.jpg']);
                imsz = size(im);
            catch ne_eo;
                imsz = [1, 1, 3];
            end
        end
        csz.(msm) = imsz(1:2);
    end
    imsz = csz.(msm);
    row = 1 + imsz(1) * (1 - ms(3));
    col = 1 + imsz(2) * ms(2);
    
    % determine pixel spacing
    txpix = txcrd * [imsz(2), 0; 0, imsz(1)];
    txdst = sqrt(sum((txpix([1,1], :) - txpix(2:3, :)) .^ 2, 2));
    ppmm = mean(txdst ./ trdst);
    
    % define coordinate
    row = floor(row);
    col = floor(col);
    boxhw = ceil(opts.boxwidth * 0.5 * icc.majorAxisMM * ppmm);
    x1 = max(1, col - boxhw);
    x2 = min(imsz(2), col + boxhw);
    y1 = max(1, row - boxhw);
    y2 = min(imsz(1), row + boxhw);
    
    % fill row
    lr = l(c, :);
    lr(1:nf) = {[msm '.jpg'], sel(c), imsz(2), imsz(1), x1, x2, y1, y2};
    for fc = 1:no
        lr{nf + fc} = icc.(opts.output{fc});
    end
    l(c, :) = format_values(lr);
end

% remove empty?
if opts.delempty
    l(cellfun(@isempty, l(:,1)), :) = [];
end

% add header
l = [[fo', opts.output']; l];


% sub function for value formating
function r = format_values(r)
for c = 1:numel(r)
    if ischar(r{c})
        r{c} = deblank(r{c});
        continue;
    end
    r{c} = deblank(sprintf('%g', r{c}));
end
