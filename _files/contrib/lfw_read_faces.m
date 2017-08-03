function [faces, imc, imi, imm, mimc, sdimc, opts] = lfw_read_faces(opts)
% lfw_read_faces  - read Labeled Faces in the Wild database
%
% FORMAT:       [faces, imc, imi, imm, mimc, sdimc, opts] = lfw_read_faces([opts])
%
% Input fields:
%
%       opts        optional settings
%        .folder    folder where LFW database resides (default: pwd)
%        .mask      mask images, either 250x250 logical or double or
%                   1x3 value, relative SD threshold, FWHM, and FWHM-thresh
%        .maskrgb   1x3 value placed into masked voxels (default: [0, 0, 0])
%        .random    1x1 double, read random selection of N faces
%        .reborder  replace black border with adjacent pixels (default: true)
%        .unique    read only one face per labeled person (default: true)
%
% Output fields:
%
%       faces       list of face filenames
%       imc         250x250x3xN uint8 images content
%       imi         image index (into sorted 13,233 database list)
%       imm         mean of images (after working on them)
%       mimc        mean of images (prior to any alterations)
%       sdimc       SD of images (prior to any alterations)
%       opts        opts output (for mask image, etc.)
%
% Note: this function requires the deep-funneled aligned LFW database with
%       13,233 faces (5,749 unique).
%
% Note: given the requirements for memory (~2.5GB for all faces in uint8)
%       this function should only be used (for the full dataset) with
%       at least 12GB of free memory (all images, flexi-mask option used)
%       and at least 5GB of free memory (unique images, flex-mask option)
%
% Example:
%
% [faces, imc, imi, imm, mimc, sdimc] = lfw_read_faces(struct( ...
%     'mask', [0.85, 10, 12]))

% requires neuroelf
n = neuroelf;

% options
if nargin < 1 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'folder') || ...
   ~ischar(opts.folder) || ...
    isempty(opts.folder) || ...
    exist(opts.folder(:)', 'dir') ~= 7
    opts.folder = pwd;
end
if ~isfield(opts, 'mask') || ...
   ((~isa(opts.mask, 'double') || ...
     (~isequal(size(opts.mask), [1, 3]) && ...
      ~isequal(size(opts.mask), [250, 250]))) && ...
    (~islogical(opts.mask) || ...
     ~isequal(size(opts.mask), [250, 250])))
    opts.mask = [];
elseif numel(opts.mask) == 3 && ...
   (any(isinf(opts.mask) | isnan(opts.mask) | opts.mask <= 0, 2) || ...
    opts.mask(1) >= 1 || ...
    opts.mask(2) < 0.5 || ...
    opts.mask(3) > 1000)
    opts.mask = [];
end
if ~isfield(opts, 'maskrgb') || ...
   ~isa(opts.maskrgb, 'double') || ...
    numel(opts.maskrgb) ~= 3 || ...
    any(isinf(opts.maskrgb) | isnan(opts.maskrgb) | opts.maskrgb < 0 | opts.maskrgb > 255)
    opts.maskrgb = [0, 0, 0];
else
    opts.maskrgb = round(opts.maskrgb(:))';
end

% there are a total of 13,233 faces in the database
if ~isfield(opts, 'random') || ...
   ~isa(opts.random, 'double') || ...
    numel(opts.random) ~= 1 || ...
    isinf(opts.random) || ...
    isnan(opts.random) || ...
    opts.random < 1 || ...
    opts.random > 13233
    opts.random = [];
else
    opts.random = round(opts.random);
end
if ~isfield(opts, 'reborder') || ...
   ~islogical(opts.reborder) || ...
    numel(opts.reborder) ~= 1
    opts.reborder = true;
end
if ~isfield(opts, 'unique') || ...
   ~islogical(opts.unique) || ...
    numel(opts.unique) ~= 1
    opts.unique = true;
    
    % there are only 5,749 unique faces in the database
    if ~isempty(opts.random) && ...
        opts.random > 5749
        opts.random = [];
    end
end

% find faces
faces = n.findfiles(opts.folder, '*.jpg');
nfaces = numel(faces);

% get folders (IDs) of faces
faceids = n.mfileparts(faces);
nfaceids = numel(unique(faceids));

% complete database?
if nfaces ~= 13233 || ...
    nfaceids ~= 5749
    warning( ...
        'LFW:IncompleteDatabase', ...
        'Incomplete/corrupt faces LFW database (%d faces, %d unique).', ...
        nfaces, nfaceids ...
    );
    if ~isempty(opts.random)
        if opts.unique && ...
            opts.random > nfaceids
            opts.random = [];
        elseif opts.random > nfaces
            opts.random = [];
        end
    end
end

% create array to hold data
imc = uint8(zeros(250, 250, 3));

% non-random
if isempty(opts.random)
    
    % unique
    if opts.unique
        
        % get unique list
        ridx = randperm(nfaces);
        [ufaces, uid] = unique(faceids(ridx), 'stable');
        imi = sort(ridx(uid));
        
    % all faces
    else
        
        % take indices
        imi = (1:nfaces)';
    end
    
% random selection (and order)
else
    
    % make random selection
    ridx = randperm(nfaces);
    
    % if unique faces
    if opts.unique
        
        % get unique list of face names
        [ufaces, uid] = unique(faceids(ridx), 'stable');
        imi = ridx(uid(1:opts.random));
        
    % non-unique faces
    else
        
        % just take first N
        imi = ridx(1:opts.random);
    end
end

% extend array
nimi = numel(imi);
imc(1, 1, 1, nimi) = 0;

% read faces
dc = 0.019;
di = 0.2;
fprintf('reading face images: ');
for fc = 1:nimi
    try
        sim = imread(faces{imi(fc)});
        if size(sim, 3) == 1
            sim = sim(:, :, [1, 1, 1]);
        end
        imc(:, :, :, fc) = sim;
        if fc > (dc * nimi)
            fprintf('.');
            pause(0.01);
            dc = dc + 0.02;
            if dc > di
                fprintf(' %d%% ', round(100 * di));
                di = di + 0.2;
            end
        end
    catch ne_eo;
        warning( ...
            'LFW:ImageReadError', ...
            ne_eo.message ...
        );
    end
end
fprintf('\n');

% compute mean, SD (over non-0 pixels)
nimc = zeros(250, 250);
mimc = zeros(250, 250, 3);
sdimc = zeros(250, 250, 3);
dc = 0.019;
di = 0.2;
fprintf('computing mean & SD: ');
for fc = 1:nimi
    bimc = double(any(imc(:, :, :, fc) > 4, 3) | all(imc(:, :, :, fc) > 2, 3));
    nimc = nimc + bimc;
    bimc = repmat(bimc, [1, 1, 3]) .* double(imc(:, :, :, fc));
    mimc = mimc + bimc;
    sdimc = sdimc + bimc .* bimc;
    if fc > (dc * nimi)
        fprintf('.');
        pause(0.01);
        dc = dc + 0.02;
        if dc > di
            fprintf(' %d%% ', round(100 * di));
            di = di + 0.2;
        end
    end
end
sdimc = (sdimc - (mimc .* mimc) ./ repmat(nimc, [1, 1, 3])) ./ (repmat(nimc - 1, [1, 1, 3]));
sdimc(isinf(sdimc) | isnan(sdimc) | sdimc < 0) = 0;
mimc = mimc ./ repmat(nimc, [1, 1, 3]);
msdimc = max(sdimc(:));
sdimc(sdimc == 0) = msdimc;
sdimc = sqrt(sdimc);
fprintf('\n');

% 1x3 mask
if numel(opts.mask) == 3
    
    % compute mean of SD along RGB dim
    msdimc = mean(sdimc, 3);
    
    % make binary mask
    immsk = double(msdimc <= (opts.mask(1) * max(msdimc(:))));
    
    % smooth mask
    smk = n.smoothkern(opts.mask(2));
    smk(smk < 0.0001) = [];
    immsk = n.flexinterpn(immsk, [Inf, Inf; ones(2, 2); 250, 250], smk, 1, 0);
    
    % threshold mask
    opts.mask = min(1, immsk ./ (opts.mask(3) .* (max(smk) ^ 2)));
end

% apply mask
if ~isempty(opts.mask)
    
    % binary mask
    if islogical(opts.mask)
        opts.imask = ~opts.mask;
        
    % double mask
    else
        opts.mask = n.limitrangec(opts.mask, 0, 1, 0);
        opts.imask = 1 - opts.mask;
    end
    
    % iterate
    dc = 0.019;
    di = 0.2;
    fprintf('masking face images: ');
    for fc = 1:nimi
        
        % binary mask
        if islogical(opts.mask)
            
            % RGB planes
            for pc = 1:3
                p = imc(:, :, pc, fc);
                p(opts.imask) = opts.maskrgb(pc);
                imc(:, :, pc, fc) = p;
            end
            
        % double mask
        else
            
            % RGB planes
            for pc = 1:3
                imc(:, :, pc, fc) = uint8(round( ...
                    double(imc(:, :, pc, fc)) .* opts.mask + ...
                    opts.maskrgb(pc) .* opts.imask));
            end
        end
        if fc > (dc * nimi)
            fprintf('.');
            pause(0.01);
            dc = dc + 0.02;
            if dc > di
                fprintf(' %d%% ', round(100 * di));
                di = di + 0.2;
            end
        end
    end
    fprintf('\n');
end

% compute mean
imm = mean(imc, 4);
