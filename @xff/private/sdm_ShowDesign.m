function f = sdm_ShowDesign(xo, opts)
% SDM::ShowDesign  - show design matrix graphically
%
% FORMAT:       [f = ] sdm.ShowDesign([opts])
%
% Input fields:
%
%       opts        optional settings
%        .confounds flag, also include confounds in display, default: false
%        .type      either of {'image'}, 'ortho', 'plot'
%
% Output fields:
%
%       f           figure handle of plot (Matlab handle!)
%
% Using: orthimage.

% Version:  v1.1
% Build:    16021413
% Date:     Feb-14 2016, 1:02 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2014, 2016, Jochen Weber
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

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'sdm')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'confounds') || ~islogical(opts.confounds) || numel(opts.confounds) ~= 1
    opts.confounds = false;
end
if ~isfield(opts, 'type') || ~ischar(opts.type) || ...
   ~any(strcmpi(opts.type(:)', {'i', 'image', 'o', 'ortho', 'p', 'plot'}))
    opts.type = 'i';
else
    opts.type = lower(opts.type(1));
end

% get design matrix right
bc = xo.C;
fc = bc.FirstConfoundPredictor;
sdm = bc.SDMMatrix;
sdmf = xo.F;
if isempty(sdmf)
    sdmf = 'interactive.sdm';
end
[sdmp, sdmf] = fileparts(sdmf);

% what to plot
if ~opts.confounds
    if fc >= size(sdm, 2)
        sdm(:, fc:end) = [];
    end
end

% create new figure and axes
f = figure;
x = axes('Parent', f);

% how to plot
switch (opts.type)

    % image
    case 'i'

        % create gray-scale version
        msdm = max(abs(sdm(:)));
        gsdm = floor(min(64 - sqrt(eps), max(0, 32 + (32 / msdm) .* sdm)));
        image(gsdm, 'Parent', x);
        colormap(x, 'gray');
        dtype = 'Design matrix (image)';

    % orthogonality
    case 'o'

        % plot orthogonality image
        sdmoi = floor((64 - sqrt(eps)) .* ne_methods.orthimage(sdm));
        image(sdmoi, 'Parent', x);
        colormap(x, 'gray');
        dtype = 'Design orthogonality';

    % time-course plot
    case 'p'

        % plot with distance of 2 between regressors
        plot(sdm + 2 .* ones(size(sdm, 1), 1) * ((size(sdm, 2) - 1):-1:0));
        dtype = 'Design matrix (time courses)';
        legend(bc.PredictorNames{1:size(sdm, 2)});
end

% set figure name
set(f, 'NumberTitle', 'off', 'Name', sprintf('%s: %s, %d-x-%d', dtype, sdmf, size(sdm)));
