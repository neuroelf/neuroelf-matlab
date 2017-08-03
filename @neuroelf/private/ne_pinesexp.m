% PUBLIC FUNCTION ne_pinesexp: run PINES expression analysis
function varargout = ne_pinesexp(varargin)

% Version:  v1.1
% Build:    17062813
% Date:     Jun-28 2017, 1:30 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2016, 2017, Jochen Weber
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
cc = ne_gcfg.fcfg;
ch = ne_gcfg.h;
mFig = ch.MainFig;
mfp = mFig.Pointer;
pbar = ch.Progress;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% only allow one instance
if any(strcmp(ne_gcfg.c.blockcb, 'pinesexp'))
    return;
end
ne_gcfg.c.blockcb{end+1} = 'pinesexp';

% requires current statsvar to be GLM
glm = cc.StatsVar;
if numel(glm) ~= 1 || ~isxff(glm, 'glm') || glm.ProjectType ~= 1
    glm = [];
    if isempty(glm) && isfield(cc, 'CM') && isstruct(cc.CM) && isfield(cc.CM, 'glm') && ...
        numel(cc.CM.glm) == 1 && isxff(cc.CM.glm, 'glm') && cc.CM.glm.ProjectType == 1
        glm = cc.CM.glm;
    end
end
if numel(glm) ~= 1 || ~isxff(glm, 'glm') || glm.ProjectType ~= 1
    uiwait(warndlg('PINES expression analysis requires a GLM to be loaded.', 'NeuroElf - info', 'modal'));
    ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'pinesexp')) = [];
    return;
end
bbox = glm.BoundingBox;
rtv = glm.RunTimeVars;
cons = rtv.Contrasts;
if isempty(cons)
    cons = cell(0, 2);
end
ncons = size(cons, 1);
subs = glm.Subjects;
nsubs = numel(glm.Subjects);
spred = glm.SubjectPredictors;
if strcmpi(spred{end}, 'constant')
    spred(end) = [];
end
stpred = ~cellfun('isempty', regexpi(spred, '_T\d+$'));
stpreds = (sum(stpred) >= 30 && (nargin < 3 || ~isa(varargin{3}, 'double') || isempty(varargin{3}) || ...
    any(isinf(varargin{3}(:)) | isnan(varargin{3}(:)) | varargin{3}(:) < 1)));
nspred = numel(spred);
tspred = nspred;
if stpreds
    spredi = spred;
    for pc = 1:nspred
        if ~stpred(pc)
            spredi{pc} = pc;
        elseif ~isempty(regexpi(spred{pc}, '_T0+1$'))
            spredi{pc} = lsqueeze(find(~cellfun('isempty', ...
                regexp(spred, ['^' regexprep(spred{pc}, '_T0+1$', '_T\\d+$')]))));
        else
            spredi{pc} = zeros(0, 1);
        end
    end
    kspred = find(~cellfun('isempty', spredi));
    spred = regexprep(spred(kspred), '_T0+1$', '_T###');
    nspred = numel(spred);
end
ssel = rtv.SubSels;
if ~isempty(ssel)
    ssel = multimatch(ssel{1, 2}, subs);
else
    ssel = lsqueeze(1:nsubs)';
end
cons(ncons+1:ncons+numel(spred), 1) = spred;
for pc = 1:nspred
    cvals = zeros(tspred, 1);
    if stpreds
        cvals(spredi{kspred(pc)}) = 1;
    else
        cvals(pc) = 1;
    end
    cons{ncons+pc, 2} = cvals;
end

% what contrast
if nargin < 3 || ~isa(varargin{3}, 'double') || isempty(varargin{3}) || ...
    any(isinf(varargin{3}(:)) | isnan(varargin{3}(:)) | varargin{3}(:) < 1)
    [cidx, cok] = listdlg('ListString', cons(:, 1), 'SelectionMode', 'multiple', ...
        'ListSize', [min(640, max(320, 10 * size(char(cons(:, 1)), 2))), 12 * (size(cons, 1) + 2)], ...
        'InitialValue', [], 'Name', 'NeuroElf - selection', ...
        'PromptString', 'Please select the maps for which you wish to compute the pattern expression score...');
    if isempty(cok) || isempty(cidx) || isequal(cok, 0)
        ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'pinesexp')) = [];
        return;
    end
else
    cidx = unique(round(min(size(cons, 1), varargin{3}(:))));
end
cons = cons(cidx, :);
if stpreds
    xcons = find(cidx > ncons);
    for pc = numel(xcons):-1:1
        xconc = xcons(pc);
        scon = sum(cons{xconc, 2});
        if scon > 1
            fcon = find(cons{xconc, 2}(:));
            cname = regexprep(cons{xconc, 1}, '_T\#+$', '_T');
            cons = [cons(1:xconc-1, :); cell(scon, 2); cons(xconc+1:end, :)];
            for pcc = 1:scon
                cons{xconc+pcc-1, 1} = sprintf('%s%03d', cname, pcc);
                cvals = zeros(tspred, 1);
                cvals(fcon(pcc)) = 1;
                cons{xconc+pcc-1, 2} = cvals;
            end
        end
    end
end
ncons = size(cons, 1);

% subject selection
if nargin < 4 || ~isa(varargin{4}, 'double') || isempty(varargin{4}) || ...
    any(isinf(varargin{4}(:)) | isnan(varargin{4}(:)) | varargin{4}(:) < 1)
    [sidx, cok] = listdlg('ListString', subs, 'SelectionMode', 'multiple', ...
        'ListSize', [480, min(640, max(320, 12 * nsubs + 2))], ...
        'InitialValue', ssel, 'Name', 'NeuroElf - selection', ...
        'PromptString', 'Please select the subjects for which you wish to compute the pattern expression score...');
    if isempty(cok) || isempty(sidx) || isequal(cok, 0)
        ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'pinesexp')) = [];
        return;
    end
else
    sidx = unique(round(min(numel(subs), varargin{4}(:))));
end
ssel = subs(sidx);
nsubs = numel(sidx);

% get a few parameters
params = inputdlg( ...
    {'Computation: (c)orrelation, or (d)ot-product?', 'Apply addition mask: (y)es or (n)o?', ...
     'Sort conditions alphabetically: (yes) or (n)o?'}, 'PINES expression analysis configuration', ...
    1, {'  d', '  n', '  n'});
if numel(params) ~= 3 || ~iscell(params)
    ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'pinesexp')) = [];
    return;
end
params = ddeblank(params(:));
if any(cellfun('prodofsize', params) ~= 1)
    ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'pinesexp')) = [];
    return;
end

% sort conditions
if lower(params{3}) == 'y'
    [cons(:, 1), csort] = sort(cons(:, 1));
    cons(:, 2) = cons(csort, 2);
end

% apply additional mask
if lower(params{2}) == 'y'
    mskfile = xff('*.(hdr|msk|vmr)', 'Please select a file to mask the PINES expression analysis with...');
    if isempty(mskfile) || ~isxff(mskfile, true)
        ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'pinesexp')) = [];
        return;
    end
    mskfile.LoadTransIOData;
    if isxff(mskfile, 'hdr')
        mskfile.VoxelData = single(mskfile.VoxelData ~= 0);
    end
    mskbox = mskfile.SampleBVBox(bbox, 1, 'linear');
    mskbox = (~isinf(mskbox) & ~isnan(mskbox) & mskbox >= 0.5);
    mskfile.ClearObject;
else
    mskbox = [];
end

% computation type
usecorrel = (lower(params{1}) == 'c');
if usecorrel
    ctext = 'correlation';
else
    ctext = 'dot-product';
end

% output data
pinesval = zeros(nsubs, ncons);
dlines = cell(nsubs + 1, ncons + 1);
dlines{1} = 'subject';
dlines(1, 2:end) = cons(:, 1)';
dlines(2:end, 1) = ssel(sidx);

% pointer
mFig.Pointer = 'watch';
tprog = 1 + 2 * numel(pinesval);
cprog = ne_progress(0, 0, {true, 0, 'Sampling and masking pattern image...'});
nprog = 1 / tprog;
drawnow;

% load and sample PINES map
if nargin > 4 && ischar(varargin{5}) && strcmpi(varargin{5}(:)', 'select')
    [pinesfile, pinespath] = uigetfile({'*.hdr;*.nii;*.nii.gz;*.head;*.mgh;*.mgz;*.vmp', ...
        'All compatible pattern files (*.hdr, *.nii, *.head, *.mgh, *.vmp)'}, ...
        'Please select a compatible pattern file...');
    if isequal(pinesfile, 0) || isequal(pinespath, 0) || isempty(pinesfile)
        ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'pinesexp')) = [];
        return;
    end
    if isempty(pinespath)
        pinespath = pwd;
    end
    pinesfile = [pinespath filesep pinesfile];
    if exist(pinesfile, 'file') == 2
        varargin{5} = pinesfile;
    else
        ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'pinesexp')) = [];
        return;
    end
end
if nargin < 5 || ~ischar(varargin{5}) || isempty(varargin{5}) || exist(varargin{5}(:)', 'file') ~= 2
    pines = neuroelf_file('n', 'PINES_PCR');
    pinesfile = 'PINES image';
else
    try
        pines = [];
        pinesfile = varargin{5}(:)';
        pines = xff(pinesfile);
        if ~isxff(pines, {'hdr', 'head', 'mgh', 'vmp'})
            error('neuroelf:xff:wrongObject', 'Invalid pattern file.');
        end
    catch ne_eo;
        if isxff(pines, true)
            pines.ClearObject;
        end
        uiwait(warndlg(['Error loading pattern image: ' ne_eo.message], ...
            'NeuroElf - error', 'modal'));
        ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'pinesexp')) = [];
        return;
    end
end
pines.LoadTransIOData;
pinesbox = pines.SampleBVBox(bbox, 1, cc.imethod);
if isxff(pines, 'hdr')
    pines.VoxelData = single(pines.VoxelData ~= 0);
elseif isxff(pines, 'vmp')
    pines.Map(1).VMPData = single(pines.Map(1).VMPData ~= 0);
elseif isxff(pines, 'head')
    pines.Brick(1).Data = single(pines.Brick(1).Data ~= 0);
end
pinesmsk = pines.SampleBVBox(bbox, 1, 'linear');
pinesmsk = (~isinf(pinesmsk) & ~isnan(pinesmsk) & pinesmsk >= 0.5);
if ~isempty(mskbox)
    pinesmsk(~mskbox) = false;
end
pinesmsk = pinesmsk(:);
pines.ClearObject;
msize = numel(pinesbox);
pinesbox = pinesbox(pinesmsk);
npinesvx = numel(pinesbox);

% get contrasts and compute
for conc = 1:ncons
    cicd = (sum(cons{conc, 2}) == 1 && sum(cons{conc, 2} == 1) == 1);
    if cicd
        pbar.Progress(nprog, sprintf('Sampling condition %s...', cons{conc, 1}));
    else
        pbar.Progress(nprog, sprintf('Sampling contrast %s...', cons{conc, 1}));
    end
    nprog = nprog + nsubs / tprog;
    conmaps = glm.RFX_conmaps(cons{conc, 2}, struct('subsel', sidx));
    conmaps = reshape(conmaps, msize, nsubs);
    conmaps = conmaps(pinesmsk, :);
    conmaps(:, any(isnan(conmaps), 1)) = NaN;
    conmask = (~any(isinf(conmaps), 2) & (sum(conmaps ~= 0, 2) >= (0.5 * nsubs)));
    for subc = 1:nsubs
        nprog = nprog + 1 / tprog;
        if usecorrel
            cval = corrcoef([conmaps(conmask, subc), pinesbox(conmask)]);
            pinesval(subc, conc) = cval(2);
        else
            pinesval(subc, conc) = sum(conmaps(conmask, subc) .* pinesbox(conmask));
        end
        dlines{subc+1, conc+1} = sprintf('%.8g', pinesval(subc, conc));
    end
end
dlines = dlines';
dlines = sprintf([repmat('%s\t', 1, conc + 1) '\n'], dlines{:});

% display results
if nargin < 4
    ch.ClusterTable.String = ...
        sprintf('PINES expression analysis (%d subjects, %d contrasts, %s of %d voxels)\n\n%s', ...
        nsubs, ncons, ctext, npinesvx, dlines);
    assignin('base', 'PINES_scores', pinesval);
else
    fprintf('Expression analysis (%d subjects, %d contrasts, %s of %d voxels)\nImage: %s\n\n%s\n', ...
        nsubs, ncons, ctext, npinesvx, pinesfile, dlines);
    if nargout > 1
        varargout{2} = pinesfile;
    end
end

% re-enable next analysis run
ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'pinesexp')) = [];

% done
ne_progress(0, 0, cprog);
mFig.Pointer = mfp;
drawnow;

% return
if nargout > 0
    varargout{1} = pinesval;
end
