function images2movie(srcfiles, targetfile, opts)
% images2movie  - read images and turn into movie (AVI) file
%
% FORMAT:       images2movie(srcfiles, target [, opts])
%
% Input fields:
%
%       srcfiles    cell array with images (readable with imread)
%       target      target filename
%       opts        optional settings
%        .fps       frames-per-second (default: 30)
%        .pbar      either xprogress or xfigure:XProgress object
%        .profile   one of the supported profiles for VideoWriter objects
%        .rect      1x4 rectangle definition to crop image/s
%        .resize    1x2 resizing of (cropped) image to target size
%
% No output fields.
%
% Note: in addition to the options, all other fields beginning with a
%       capital letter will be attempted to set in the VideoWriter object
%       after it is created.

% Version:  v0.9d
% Build:    14062215
% Date:     Jun-22 2014, 3:20 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2014, Jochen Weber
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
if nargin < 2 || ...
   ~iscell(srcfiles) || ...
    numel(srcfiles) < 2 || ...
   ~ischar(targetfile) || ...
    numel(targetfile) < 5
    return;
end
targetfile = targetfile(:)';
if nargin < 3 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'fps') || ...
   ~isa(opts.fps, 'double') || ...
    numel(opts.fps) ~= 1 || ...
    isinf(opts.fps) || ...
    isnan(opts.fps) || ...
    opts.fps < 1
    opts.FrameRate = 30;
else
    opts.FrameRate = min(120, opts.fps);
end
if ~isfield(opts, 'pbar') || ...
    numel(opts.pbar) ~= 1 || ...
   ~any(strcmpi(class(opts.pbar), {'xfigure', 'xprogress'}))
    opts.pbar = [];
end
if ~isfield(opts, 'profile') || ...
   ~ischar(opts.profile) || ...
    isempty(opts.profile) || ...
   ~any(strcmp(opts.profile(:)', ...
        {'Archival', 'Motion JPEG AVI', 'Motion JPEG 2000', 'MPEG-4', ...
         'Uncompressed AVI', 'Indexed AVI', 'Grayscale AVI'}))
    opts.profile = 'MPEG-4';
else
    opts.profile = opts.profile(:)';
end
if ~isfield(opts, 'rect') || ...
   ~isa(opts.rect, 'double') || ...
    numel(opts.rect) ~= 4 || ...
    any(isinf(opts.rect(:)) | isnan(opts.rect(:)) | opts.rect(:) < 1 | opts.rect(:) ~= round(opts.rect(:)))
    opts.rect = [];
else
    opts.rect = opts.rect(:)';
    opts.rect(3:4) = max(opts.rect(1:2), opts.rect(3:4));
end
if ~isfield(opts, 'resize') || ...
   ~isa(opts.resize, 'double') || ...
    numel(opts.resize) ~= 2 || ...
    any(isinf(opts.resize) | isnan(opts.resize) | opts.resize < 1)
    opts.resize = [];
else
    opts.resize = ceil(opts.resize(:)');
end

% read first image
try
    im = imread(srcfiles{1});
catch ne_eo;
    rethrow(ne_eo);
end

% rectangle
if isempty(opts.rect)
    opts.rect = [1, 1, size(im, 2), size(im, 1)];
end

% test xprogress
nims = numel(srcfiles);
if nims >= opts.FrameRate
    pst = 1 / 86400;
    pbarn = now + pst;
    try
        if isempty(opts.pbar)
            pbar = xprogress;
            xprogress(pbar, 'setposition', [80, 200, 640, 36]);
            xprogress(pbar, 'settitle', 'Converting images to movie...');
            xprogress(pbar, 0, sprintf('0/%d images written.', nims), 'visible', 0, 1);
        else
            pbar = opts.pbar;
            pbarvis = pbar.Visible;
            pbar.Progress(0, sprintf('Creating movie from %d images...', nims));
            pbar.Visible = 'on';
        end
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        pbar = [];
        pbarn = Inf;
    end
else
    pbar = [];
    pbarn = Inf;
end

% big try loop
try
    
    % create video stream
    m = VideoWriter(targetfile);
    
    % settings
    of = fieldnames(opts);
    for fc = 1:numel(of)
        if ~isempty(of{fc}) && ...
            of{fc}(1) == upper(of{fc}(1)) && ...
           ~any(strcmpi(of{fc}, {'fps', 'pbar', 'profile', 'rect', 'resize'}))
       
            % without a fuss...
            try
                m.(of{fc}) = opts.(of{fc});
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
            end
        end
    end
    
    % open file
    open(m);
    
    % start putting in content
    for fc = 1:nims
        
        % read image
        im = imread(srcfiles{fc});
        
        % rectangle
        im = im(opts.rect(2):opts.rect(4), opts.rect(1):opts.rect(3), :);
        
        % grayscale -> RGB
        if size(im, 3) == 1
            im = im(:, :, [1, 1, 1]);
        end
        
        % resize
        if ~isempty(opts.resize)
            im = image_resize(im, opts.resize(2), opts.resize(1));
        end
        
        % write into video stream
        cf = struct('cdata', im, 'colormap', []);
        writeVideo(m, cf);
        
        % progress bar
        if pbarn <= now
            if isempty(opts.pbar)
                pbar.Progress(fc / nims, sprintf('%d/%d images written.', fc, nims));
            else
                pbar.Progress(fc / nims);
            end
            pbarn = now + pst;
        end
    end
    if ~isempty(pbar)
        pbar.Progress(1);
    end
    
    % close video
    close(m);
    
catch ne_eo;
    le = ne_eo;
    try
        close(m);
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
    end
    rethrow(le);
end

% progress bar
if ~isempty(pbar)
    if isempty(opts.pbar)
        closebar(pbar);
    else
        pbar.Visible = pbarvis;
    end
end
