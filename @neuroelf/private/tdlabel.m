function [tl, tcs] = tdlabel(c, opts)
% tdlabel  - get NGM label of coordinate
%
% FORMAT:       label = tdlabel(c [, opts])
%
% Input fields:
%
%       c           Cx3 coordinate
%       opts        optional fields
%        .dist      give distance if no direct hit, default: true
%        .gyri      give structure if no direct area name, default: true
%        .hemi      give hemisphere in RH/LH notation, default: true
%        .hitcoord  give coordinate hit if no direct hit, default: true
%        .maxdist   max NGM distance, default: 10mm
%
% Output fields:
%
%       label       Cx1 cell array with labels

% Version:  v0.9a
% Build:    11031817
% Date:     May-17 2010, 10:48 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, Jochen Weber
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

% check arguments
if nargin < 1 || ...
   ~isa(c, 'double') || ...
    size(c, 2) ~= 3 || ...
    isempty(c) || ...
    any(isinf(c(:)) | isnan(c(:)))
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid/missing coordinates in call.' ...
    );
end
if nargin < 2 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'dist') || ...
   ~islogical(opts.dist) || ...
    numel(opts.dist) ~= 1
    opts.dist = true;
end
if ~isfield(opts, 'gyri') || ...
   ~islogical(opts.gyri) || ...
    numel(opts.gyri) ~= 1
    opts.gyri = true;
end
if ~isfield(opts, 'hemi') || ...
   ~islogical(opts.hemi) || ...
    numel(opts.hemi) ~= 1
    opts.hemi = true;
end
if ~isfield(opts, 'hitcoord') || ...
   ~islogical(opts.hitcoord) || ...
    numel(opts.hitcoord) ~= 1
    opts.hitcoord = true;
end
if ~isfield(opts, 'maxdist') || ...
   ~isa(opts.maxdist, 'double') || ...
    numel(opts.maxdist) ~= 1 || ...
    isinf(opts.maxdist) || ...
    isnan(opts.maxdist) || ...
    opts.maxdist < 1 || ...
    opts.maxdist > 20
    opts.maxdist = 10;
end

% build output array
nc = size(c, 1);
tl = cell(nc, 1);
tcs = c;

% iterate over coordinates
for cc = 1:nc

    % call tdlocal2
    tc = c(cc, :);
    ngm = tdlocal2(-5, tc(1), tc(2), tc(3), [opts.maxdist, 1]);

    % make some checks
    if ~iscell(ngm)
        if ~isempty(ngm)
            ngm = {{tc, [0, 0, 0], 0, ngm}};
        else
            ngm = {{[], [0, 0, 0], 0, {'*', '*', '*', '*', '*'}}};
        end
    end

    % get first (and only) hit
    ngm = ngm{1};

    % extract label informations
    if ~iscell(ngm{4})
        ngm{4} = splittocell(ngm{4}, ',');
    end

    % only accept full labels
    if numel(ngm{4}) ~= 5
        ngm{4} = {'*', '*', '*', '*', '*'};
    end
    tlp = ngm{4};

    % only produce label if anything found
    if ~isempty(ngm{1}) && ...
       ~strcmp(tlp{1}, '*')

        % store coordinate
        tcs(cc, :) = ngm{1};

        % extract hemisphere
        if opts.hemi
            talhemi = splittocell(tlp{1}, ' ');
            if isempty(talhemi)
                talhemi = '';
            else
                talhemi = [talhemi{1}(1) 'H '];
            end
        else
            talhemi = '';
        end

        % get gyri information
        if ~isempty(tlp{5}) && ...
           ~strcmp(tlp{5}, '*') && ...
            opts.gyri
            if ~isempty(tlp{3}) && ...
               ~strcmp(tlp{3}, '*')
                talgyri = [' (' tlp{5} ')'];
            else
                talgyri = tlp{5};
            end
        else
            talgyri = '';
        end

        % add distance information
        if ngm{3} > 0 && ...
            opts.dist

            % add hit coordinate as well
            if opts.hitcoord
                taldist = ...
                    sprintf(' (%3d;%4d;%3d) [d=%3.1fmm]', ngm{1}, ngm{3});

            % dist only
            else
                taldist = sprintf(' [d=%3.1fmm]', ngm{3});
            end

        % no distance information but hit coord all the same
        elseif ngm{3} && ...
            opts.hitcoord
            taldist = sprintf(' (%3d;%4d;%3d)', ngm{1});

        % nothing at all
        else
            taldist = '';
        end

        % structure
        if ~isempty(tlp{3}) && ...
           ~strcmp(tlp{3}, '*')
            taltext = [talhemi tlp{3} talgyri taldist];
        else
            taltext = [talhemi tlp{2} talgyri taldist];
        end

        % put into array
        tl{cc} = taltext;
    end
end
