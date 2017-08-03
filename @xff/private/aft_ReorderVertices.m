function xo = aft_ReorderVertices(xo, d)
% AFT::ReorderVertices  - reorder data in vertices between BV/NE order
%
% FORMAT:       [obj = ] obj.ReorderVertices([direction]);
%
% Input fields:
%
%       direction   one of {'auto'}, 'bv2ne', 'ne2bv'
%
% Output fields:
%
%       obj         altered object
%
% Note: this method only works for objects with 10242, 20484, 40962,
%       81924, 163842, or 327684 vertices (factor 5, 6, or 7 icosahedrons
%       and pairs of those
%
% TYPES: FSMF, GLM, MTC, SMP

% Version:  v1.1
% Build:    16031617
% Date:     Mar-16 2016, 5:27 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2014, 2015, 2016, Jochen Weber
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

% persistent data
persistent bvsph;
if isempty(bvsph)
    try
        bvsph = load([neuroelf_path('srf') filesep 'bvsph.mat']);
        srf = ne_methods.spheresrf(140, 5);
        [snei, sneim] = srf_NeighborsNDegree(srf, 1, struct('notself', true));
        delete(srf);
        sneim(1:12, :) = 1.2 .* sneim(1:12, :);
        bvsph.f5nei = sneim;
        srf = ne_methods.spheresrf(140, 6);
        [snei, sneim] = srf_NeighborsNDegree(srf, 1, struct('notself', true));
        delete(srf);
        sneim(1:12, :) = 1.2 .* sneim(1:12, :);
        bvsph.f6nei = sneim;
        srf = ne_methods.spheresrf(140, 7);
        [snei, sneim] = srf_NeighborsNDegree(srf, 1, struct('notself', true));
        delete(srf);
        sneim(1:12, :) = 1.2 .* sneim(1:12, :);
        bvsph.f7nei = sneim;
    catch xfferror
        bvsph = [];
        rethrow(xfferror);
    end
end

% check arguments
if numel(xo) ~= 1 || ~xffisobject(xo, true, {'fsmf', 'glm', 'mtc', 'smp'})
    error('neuroelf:xff:badArgument', 'Bad object or argument in call.');
end
bc = xo.C;
if strcmpi(xo.S.Extensions{1}, 'glm') && bc.ProjectType ~= 2
    error('neuroelf:xff:badArgument', 'Only valid for MTC-based GLMs.');
end
if ~any(bc.NrOfVertices == [10242, 20484, 40962, 81924, 163842, 327684])
    error('neuroelf:xff:badArgument', 'Invalid number of vertices.');
end
if mod(bc.NrOfVertices, 10) == 2
    tobv = bvsph.bvsph(1:bc.NrOfVertices);
    if bc.NrOfVertices == 10242
        wnei = bvsph.f5nei;
    elseif bc.NrOfVertices == 40962
        wnei = bvsph.f6nei;
    else
        wnei = bvsph.f7nei;
    end
else
    tobv = bvsph.bvsph(1:round(0.5 * bc.NrOfVertices));
    tobv(end+1:end+numel(tobv)) = tobv + numel(tobv);
    if bc.NrOfVertices == 20484
        wnei = bvsph.f5nei;
    elseif bc.NrOfVertices == 81924
        wnei = bvsph.f6nei;
    else
        wnei = bvsph.f7nei;
    end
    wnein = spalloc(size(wnei, 1), size(wnei, 2), 0);
    wnei = [[wnei, wnein]; [wnein, wnei]];
end
tone = zeros(size(tobv));
tone(tobv) = 1:numel(tobv);
if nargin < 2 || ~ischar(d) || isempty(d) || ~any(lower(d(1)) == 'abn')
    d = 'a';
else
    d = lower(d(1));
end
td = d;

% depending on object type
type = lower(xo.S.Extensions{1});
switch (type)

    % GLM
    case 'glm'

        % auto-detect
        if d == 'a'

            % RFX
            if bc.ProjectTypeRFX > 0
                tmap = bc.GLMData.RFXGlobalMap(:);
                if numel(unique(tmap)) == 1
                    tmap = bc.GLMData.Subject(1).BetaMaps(:, 1);
                end

            % FFX
            else
                tmap = bc.GLMData.TimeCourseMean(:, 1);
            end
            tmap = double(tmap);
            adir = detectshuffle(tmap, wnei, tobv, tone);
            if adir == 'b'
                to = tone;
            elseif adir == 'n'
                to = tobv;
            elseif adir == '0'
                return;
            else
                error('neuroelf:xff:detectionFailed', 'Couldn''t detect vertex order.');
            end
        elseif d == 'b'
            to = tone;
        else
            to = tobv;
        end

        % same again
        if bc.ProjectTypeRFX > 0
            bc.GLMData.RFXGlobalMap = bc.GLMData.RFXGlobalMap(to, :);
            for sc = 1:numel(bc.GLMData.Subject)
                bc.GLMData.Subject(sc).BetaMaps = bc.GLMData.Subject(sc).BetaMaps(to, :);
            end
            if isfield(bc.RunTimeVars, 'PerSubjectGLMs') && isstruct(bc.RunTimeVars.PerSubjectGLMs) && ...
                numel(bc.RunTimeVars.PerSubjectGLMs) == numel(bc.GLMData.Subject) && ...
                isfield(bc.RunTimeVars.PerSubjectGLMs, 'SEMap')
                for sc = 1:numel(bc.GLMData.Subject)
                    bc.RunTimeVars.PerSubjectGLMs(sc).SEMap = bc.RunTimeVars.PerSubjectGLMs(sc).SEMap(to, :);
                end
            end
        else
            bc.GLMData.MultipleRegressionR = bc.GLMData.MultipleRegressionR(to, :);
            bc.GLMData.MCorrSS = bc.GLMData.MCorrSS(to, :);
            bc.GLMData.BetaMaps = bc.GLMData.BetaMaps(to, :);
            bc.GLMData.XY = bc.GLMData.XY(to, :);
            bc.GLMData.TimeCourseMean = bc.GLMData.TimeCourseMean(to, :);
            if ~isempty(bc.GLMData.ARLag)
                bc.GLMData.ARLag = bc.GLMData.ARLag(to, :);
            end
        end

    % MTC
    case 'mtc'

        % auto-detect
        if d == 'a'
        end

    % SMP
    case {'fsmf', 'smp'}

        % iterate over maps
        for mc = 1:numel(bc.Map)

            % auto-detect
            if d == 'a'
                tmap = bc.Map(mc).SMPData(:);
                d = detectshuffle(tmap, wnei, tobv, tone);
                if d == '0'
                    return;
                elseif isempty(d)
                    error('neuroelf:xff:detectionFailed', 'Couldn''t detect vertex order.');
                end
            end

            % direction
            if d == 'b'
                bc.Map(mc).SMPData = bc.Map(mc).SMPData(tone, 1);
            else
                bc.Map(mc).SMPData = bc.Map(mc).SMPData(tobv, 1);
            end
        end
end

% set content back
xo.C = bc;


% sub-function
function adir = detectshuffle(tmap, wnei, tobv, tone)

% shuffle maps
tobvmap = tmap(tobv, :);
tonemap = tmap(tone, :);

% compute correlation matrices with average of neighbors
nmap = wnei * tmap;
tobvnmap = wnei * tobvmap;
tonenmap = wnei * tonemap;
gmap = ~(isinf(nmap) | isnan(nmap) | nmap == 0 | isinf(tmap) | isnan(tmap) | tmap == 0);
ccnull = corrcoef(tmap(gmap), nmap(gmap));
gmap = ~(isinf(tobvnmap) | isnan(tobvnmap) | tobvnmap == 0 | ...
    isinf(tobvmap) | isnan(tobvmap) | tobvmap == 0);
cctobv = corrcoef(tobvmap(gmap), tobvnmap(gmap));
gmap = ~(isinf(tonenmap) | isnan(tonenmap) | tonenmap == 0 | ...
    isinf(tonemap) | isnan(tonemap) | tonemap == 0);
cctone = corrcoef(tonemap(gmap), tonenmap(gmap));
if cctobv(2) > ccnull(2) && cctobv(2) > cctone(2)
    adir = 'n';
elseif cctone(2) > ccnull(2) && cctone(2) > cctobv(2)
    adir = 'b';
elseif ccnull(2) > cctobv(2) && ccull(2) > cctone(2)
    adir = '0';
else
    adir = '';
end
