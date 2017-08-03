function xo = srf_RecoToSphere(xo, settings)
% SRF::RecoToSphere  - perform all steps to get RECO mesh to sphere
%
% FORMAT:       srf = srf.RecoToSphere([settings])
%
% Input fields:
%
%       settings    optional settings
%        .distciter number of final distortion iterations (default: 0)
%        .force     boolean flag, don't heed filename
%        .smpsmooth curvature smoothing steps (default: [5, 20, 100])
%
% Output fields:
%
%       srf         altered object
%       densm       density map
%
% Note: this method calls srf.Smooth, srf.Inflate and srf.ToSphere

% Version:  v1.1
% Build:    16021110
% Date:     Feb-11 2016, 10:25 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2011, 2014, 2016, Jochen Weber
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

% global settings from config
global xffsngl;

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'srf')
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
if nargin < 2 || ~isstruct(settings) || numel(settings) ~= 1
    settings = struct;
end
if ~isfield(settings, 'distciter') || numel(settings.distciter) ~= 1 || ...
   ~isa(settings.distciter, 'double') || isinf(settings.distciter) || ...
    isnan(settings.distciter)
    try
        settings.distciter = xffsngl.CONF.settings.Morphing.ToSphere.DistCIterations;
    catch xfferror
        neuroelf_lasterr(xfferror);
        settings.distciter = 10000;
    end
else
    settings.distciter = floor(min(50000, max(0, settings.distciter)));
end
if ~isfield(settings, 'force') || numel(settings.force) ~= 1 || ~islogical(settings.force)
    settings.force = false;
end
if ~isfield(settings, 'smpsmooth') || ~isa(settings.smpsmooth, 'double') || numel(settings.smpsmooth) ~= 3 || ...
    any(isinf(settings.smpsmooth) | isnan(settings.smpsmooth) | settings.smpsmooth < 1 | settings.smpsmooth > 250)
    try
        settings.smpsmooth = xffsngl.CONF.settings.Morphing.Curvature.SmoothingSteps;
    catch xfferror
        neuroelf_lasterr(xfferror);
        settings.smpsmooth = [5, 20, 100];
    end
end

% check filename
if (numel(xo.F) < 9 || ~strcmpi(xo.F(end-8:end), '_reco.srf')) && ~settings.force
    error('neuroelf:xff:invalidFilename', 'Filename is not *_RECO.srf. Processing not forced.');
end
[sf{1:3}] = fileparts(xo.F);

% perform smoothing
[xo, d1] = srf_Smooth(xo);
aft_SaveAs(xo, [sf{1} '/' sf{2} 'SM.srf']);
aft_SaveAs(d1, [sf{1} '/' sf{2} 'SM_DENSITY.smp']);
delete(d1);

% create curvature maps
csmp = srf_CurvatureMap(xo);
for sc = 1:numel(settings.smpsmooth)
    smp_Smooth(csmp, xo, settings.smpsmooth(sc), 1);
end
aft_SaveAs(csmp, [sf{1} '/' sf{2} 'SM_CURVATURE.smp']);
delete(csmp);

% perform inflation
[xo, d2] = srf_Inflate(xo);
aft_SaveAs(xo, [sf{1} '/' sf{2} 'SM_INFL.srf']);
aft_SaveAs(d2, [sf{1} '/' sf{2} 'SM_INFL_DENSITY.smp']);
delete(d2);

% perform to-sphere morphing
[xo, d3] = srf_ToSphere(xo);
aft_SaveAs(xo, [sf{1} '/' sf{2} 'SM_INFL_SPHERE.srf']);
aft_SaveAs(d3, [sf{1} '/' sf{2} 'SM_INFL_SPHERE_DENSITY.smp']);
delete(d3);

% perform optional DC
if settings.distciter > 0
    [xo, d4] = srf_ToSphere(xo, settings.distciter, 0.2, 0.002);
    aft_SaveAs(xo, [sf{1} '/' sf{2} 'SM_INFL_SPHERE_DC.srf']);
    aft_SaveAs(d4, [sf{1} '/' sf{2} 'SM_INFL_SPHERE_DC_DENSITY.smp']);
    delete(d4);
end
