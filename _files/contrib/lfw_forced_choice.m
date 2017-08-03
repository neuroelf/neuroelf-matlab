function [c, im1, im2] = lfw_forced_choice(n, opts)
% lfw_forced_choice  - display series of face pairs to make a forced choice
%
% FORMAT:       [c, im1, im2] = lfw_forced_choice(n [, opts])
%
% Input fields:
%
%       n           make n number of forced choices (e.g. 100)
%       opts        optional settings (passed on to lfw_read_faces);
%                   alternatively, a 250x250x3xF uint8 images array
%                   if a struct is given, faces can also be given via the
%        .faces     250x250x3xF faces field
%        .mfaces    1x1 double, if set > 1 show two means of mfaces images
%        .rtweight  weight composite images by 1 / (1 + log(1 + RT)), true
%
% Output fields:
%
%       c           Nx4 array, columns: left(1), right(2), choice, RT
%
% Note: if mfaces is > 1, c will be a Nx4 cell array with indices in cells.

% check arguments
if nargin < 1 || ...
   ~isa(n, 'double') || ...
    numel(n) ~= 1 || ...
    isinf(n) || ...
    isnan(n) || ...
    n < 1 || ...
    n > 100000
    error( ...
        'LFW:BadArgument', ...
        'Invalid n.' ...
    );
else
    n = round(n);
end
if nargin < 2
    opts = struct;
end
if isa(opts, 'uint8') && ...
    size(opts, 1) == 250 && ...
    size(opts, 2) == 250 && ...
    size(opts, 3) == 3 && ...
    size(opts, 4) > 2
    opts = struct('faces', opts);
end
if numel(opts) ~= 1 || ...
   ~isfield(opts, 'faces')
    try
        opts.faces = lfw_read_faces(opts);
    catch ne_eo;
        rethrow(ne_eo);
    end
end
nfaces = size(opts.faces, 4);
if ~isfield(opts, 'mfaces') || ...
   ~isa(opts.mfaces, 'double') || ...
    numel(opts.mfaces) ~= 1 || ...
    isinf(opts.mfaces) || ...
    isnan(opts.mfaces) || ...
    opts.mfaces < 1
    opts.mfaces = 1;
else
    opts.mfaces = min(nfaces - 1, round(opts.mfaces));
end
if ~isfield(opts, 'rtweight') || ...
   ~islogical(opts.rtweight) || ...
    numel(opts.rtweight) ~= 1
    opts.rtweight = true;
end

% prepare choice array
if opts.mfaces > 1
    c = cell(n, 4);
else
    c = zeros(n, 4);
end

% initialize figure
r = get(0);
f = figure;
set(f, 'Units', 'pixels', 'Position', ...
    round([0.5 * r.ScreenSize(3) - 320, r.ScreenSize(4) - 400, 640, 320]));
ax1 = subplot(1, 2, 1);
ih1 = image(uint8(zeros(250, 250, 3)));
title(ax1, 'Face (1)');
set(ax1, 'Units', 'pixels', 'Position', [45, 35, 250, 250], 'XTick', [], 'YTick', []);
ax2 = subplot(1, 2, 2);
ih2 = image(uint8(zeros(250, 250, 3)));
title(ax2, 'Face (2)');
set(ax2, 'Units', 'pixels', 'Position', [345, 35, 250, 250], 'XTick', [], 'YTick', []);

% prepare composite images
im1 = zeros(250, 250, 3);
im2 = zeros(250, 250, 3);
mw = 0;
rtw = 1;

% iterate over choices to be made
for nc = 1:n
    
    % randomize two orders
    i1 = randperm(nfaces);
    i2 = randperm(nfaces);
    
    % simple selection
    if opts.mfaces < 2
        i1 = i1(1);
        if i2(1) ~= i1
            i2 = i2(1);
        else
            i2 = i2(2);
        end
        
        % select
        f1 = opts.faces(:, :, :, i1);
        f2 = opts.faces(:, :, :, i2);
        
    % means over selection
    else
        i1 = i1(1:opts.mfaces);
        i2 = i2(1:opts.mfaces);
        f1 = uint8(round(mean(opts.faces(:, :, :, i1), 4)));
        f2 = uint8(round(mean(opts.faces(:, :, :, i2), 4)));
    end
    
    % show faces
    set(ih1, 'CData', f1);
    set(ih2, 'CData', f2);
    
    % time now
    nt = now;
    
    % wait for key press (1) or (2)
    resp = input('Face (1) or (2): ', 's');
    if isempty(resp)
        break;
    end
    try
        resp = str2double(resp);
    catch
        resp = NaN;
    end
    
    % store
    rt = 86400 * (now - nt);
    if opts.mfaces > 1
        c(nc, :) = {i1, i2, resp, rt};
    else
        c(nc, :) = [i1, i2, resp, rt];
    end
    
    % add to composite
    if opts.rtweight
        rtw = 1 / (1 + log(1 + rt));
    end
    if ~isnan(resp) && ...
        any(resp == [1, 2])
        if resp == 1
            im1 = im1 + rtw .* double(f1);
            im2 = im2 + rtw .* double(f2);
        else
            im1 = im1 + rtw .* double(f2);
            im2 = im2 + rtw .* double(f1);
        end
        mw = mw + rtw;
    end
end

% compute composite
im1 = uint8(round(im1 ./ mw));
im2 = uint8(round(im2 ./ mw));

% delete figure;
set(ih1, 'CData', im1);
set(ih2, 'CData', im2);
