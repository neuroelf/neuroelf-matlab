function srfalignment(srfs, opts)
% srfalignment  - align surface according to their folding patterns
%
% FORMAT:       srfalignment(srfs [, opts])
%
% Input fields:
%
%       srfs        Sx1 cell array with SRF filenames
%       opts        optional settings
%        .curvsmps  Sx1 cell array with curvature SMP filenames,
%                   needed if SRFs are already spherical (curvature)
%        .icomesh   boolean flag, make icosahedron SRFs (default: true)
%        .nrvert    icosahedron number of vertices (default: 40962)
%        .spherize  boolean flag, spherize SRFs first (default: true)
%        .tempsmp   template curvature map (SMP object)
%        .tempsrf   template mesh (SRF object, needed if not icosahedron)
%
% No output fields.
%
% Note: if the input is RECO meshes (spherize := true), an average mesh
%       will be automatically created in the end

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:30 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2016, Jochen Weber
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
if nargin < 1 || ...
   ~iscell(srfs) || ...
    isempty(srfs) || ...
    numel(srfs) ~= max(size(srfs))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
srfs = srfs(:);
nsrf = numel(srfs);
if nargin < 2 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'curvsmps') || ...
   ~iscell(opts.curvsmps) || ...
    numel(opts.curvsmps) ~= numel(srfs)
    opts.curvsmps = cell(1, nsrf);
end
if ~isfield(opts, 'icomesh') || ...
   ~islogical(opts.icomesh) || ...
    numel(opts.icomesh) ~= 1
    opts.icomesh = true;
end
if ~isfield(opts, 'nrvert') || ...
   ~isa(opts.nrvert, 'double') || ...
    numel(opts.nrvert) ~= 1 || ...
    isinf(opts.nrvert) || ...
    isnan(opts.nrvert) || ...
   ~any([10242, 40962, 163842] == opts.nrvert)
    opts.nrvert = 40962;
end
if ~isfield(opts, 'spherize') || ...
   ~islogical(opts.spherize) || ...
    numel(opts.spherize) ~= 1
    opts.spherize = true;
end
if ~isfield(opts, 'tempsmp') || ...
    numel(opts.tempsmp) ~= 1 || ...
   ~isxff(opts.tempsmp, 'smp') || ...
    numel(opts.tempsmp.Map) ~= 4
    opts.tempsmp = [];
end
if ~isfield(opts, 'tempsrf') || ...
    numel(opts.tempsrf) ~= 1 || ...
   ~isxff(opts.tempsrf, 'srf')
    opts.tempsrf = [];
end

% check template
ts_c = [0, 0, 0];
ts_m = 0;
ts_t = [0, 0, 0];
ts_v = [0, 0, 0];
bvf = neuroelf_path('files');
if ~isempty(opts.tempsmp)
    ts_m = [opts.tempsmp.Map(1).SMPData; opts.tempsmp.Map(2).SMPData; ...
            opts.tempsmp.Map(3).SMPData; opts.tempsmp.Map(4).SMPData];
    if ~any([10242, 40962, 163842] == numel(opts.tempsmp.Map(1).SMPData))
        if ~isempty(opts.tempsrf) || ...
            numel(opts.tempsmp.Map(1).SMPData) ~= ...
                size(opts.tempsrf.VertexCoordinate, 1)
            error( ...
                'neuroelf:BadArgument', ...
                'NrOfVertices mismatch between tempsmp and tempsrf fields.' ...
            );
        end
        ts_c = opts.tempsrf.MeshCenter;
        ts_v = opts.tempsrf.VertexCoordinate;
        ts_t = opts.tempsrf.TriangleVertex;
    else
        if ~isempty(opts.tempsrf) && ...
            size(ts_m, 2) ~= size(opts.tempsrf.VertexCoordinate, 1)
            error( ...
                'neuroelf:BadArgument', ...
                'NrOfVertices mismatch between tempsmp and tempsrf fields.' ...
            );
        end
        if isempty(opts.tempsrf)
            try
                ico = [];
                switch (opts.nrvert)
                    case {10242}
                        ico = xff([bvf '/srf/sph20.srf']);
                    case {40962}
                        ico = xff([bvf '/srf/sph80.srf']);
                    case {163842}
                        ico = xff([bvf '/srf/sph320.srf']);
                end
                ts_c = ico.MeshCenter;
                ts_v = ico.VertexCoordinate;
                ts_t = ico.TriangleVertex;
                ico.ClearObject;
            catch ne_eo;
                if ~isempty(ico)
                    ico.ClearObject;
                end
                error( ...
                    'neuroelf:FileError', ...
                    'Error reading SPH template mesh: %s.', ...
                    ne_eo.message ...
                );
            end
        else
            ts_c = opts.tempsrf.MeshCenter;
            ts_v = opts.tempsrf.VertexCoordinate;
            ts_t = opts.tempsrf.TriangleVertex;
        end
    end
end
smps = opts.curvsmps;

% check that all files have the correct extension and exist
numvt = zeros(2, nsrf);
for sc = 1:nsrf
    if ~ischar(srfs{sc}) || ...
        numel(srfs{sc}) < 5 || ...
        numel(srfs{sc}) ~= size(srfs{sc}, 2) || ...
        exist(srfs{sc}, 'file') ~= 2 || ...
       ~strcmpi(srfs{sc}(end-3:end), '.srf')
        error( ...
            'neuroelf:BadArgument', ...
            'Invalid SRF filename in cell %d.', ...
            sc ...
        );
    end
    try
        sh = xff(srfs{sc}, 'h');
        numvt(:, sc) = [sh.NrOfVertices; sh.NrOfTriangles];
    catch ne_eo;
        error( ...
            'neuroelf:BadArgument', ...
            'Error reading SRF ''%s'' (%s).', ...
            srfs{sc}, ne_eo.message ...
        );
    end
    if numvt(2, sc) ~= (2 * (numvt(1, sc) - 2))
        error( ...
            'neuroelf:BadArgument', ...
            'SRF ''%s'' has bad euler characteristic.', ...
            srfs{sc} ...
        );
    end
    if ~opts.spherize
        if ~ischar(smps{sc}) || ...
            numel(smps{sc}) < 5 || ...
            numel(smps{sc}) ~= size(smps{sc}, 2) || ...
            exist(smps{sc}, 'file') ~= 2 || ...
           ~strcmpi(smps{sc}(end-3:end), '.smp')
            error( ...
                'neuroelf:BadArgument', ...
                'Invalid SMP filename in cell %d.', ...
                sc ...
            );
        end
        try
            sh = xff(smps{sc}, 'h');
            if sh.NrOfVertices ~= numvt(1, sc) || ...
                sh.NrOfMaps ~= 4
                error('Invalid SMP content.');
            end
        catch ne_eo;
            error( ...
                'neuroelf:BadArgument', ...
                'Invalid SMP file content in ''%s'' (%s).', ...
                smps{sc}, ne_eo.message ...
            );
        end
    end
end

% check templates
if ~isempty(opts.tempsrf)
end

% try spherize SRFs
osrf = cell(1, nsrf);
if opts.spherize
    for sc = 1:nsrf
        try
            srf = [];
            srf = xff(srfs{sc});
            srff = srf.FilenameOnDisk;
            if numel(srff) > 10 && ...
                strcmpi(srff(end-9:end), 'recosm.srf')
                smps{sc} = [srff(1:end-4) '_CURVATURE.smp'];
                osrf{sc} = srff;
                srf.RecoSMToSphere;
            else
                smps{sc} = [srff(1:end-4) 'SM_CURVATURE.smp'];
                osrf{sc} = [srff(1:end-4) 'SM.srf'];
                srf.RecoToSphere;
            end
            srfs{sc} = srf.FilenameOnDisk;
        catch ne_eo;
            if ~isempty(srf)
                srf.ClearObject;
            end
            error( ...
                'neuroelf:xffError', ...
                'Error inflating/spherizing mesh (%s).', ...
                ne_eo.message ...
            );
        end
    end
end

% icosahedron-map meshes
tsms = cell(1, nsrf);
if opts.icomesh
    try
        ico = [];
        switch (opts.nrvert)
            case {10242}
                ico = xff([bvf '/srf/sph20.srf']);
            case {40962}
                ico = xff([bvf '/srf/sph80.srf']);
            case {163842}
                ico = xff([bvf '/srf/sph320.srf']);
        end
        icoc = ico.MeshCenter;
        icov = ico.VertexCoordinate;
        icov = normvecs(icov - icoc(ones(1, size(icov, 1)), :));
        spm = xff('new:smp');
        spm.NrOfMaps = 4;
        spm.NrOfVertices = size(icov, 1);
    catch ne_eo;
        if ~isempty(ico)
            ico.ClearObject;
        end
        error( ...
            'neuroelf:FileError', ...
            'Error reading SPH template mesh: %s.', ...
            ne_eo.message ...
        );
    end
    for sc = 1:nsrf
        srff = srfs{sc};
        srf = xff(srff);
        smp = xff(smps{sc});
        srfc = srf.MeshCenter;
        srfv = srf.VertexCoordinate;
        srfv = normvecs(srfv - srfc(ones(1, size(srfv, 1)), :));
        srft = srf.TriangleVertex;
        [srfn, badsrfn, srftr] = mesh_trianglestoneighbors(size(srfv, 1), srft);
        [t, vl] = mesh_trimapmesh(icov, srfv, srft, srftr);
        if any(t == 0)
            error( ...
                'neuroelf:InternalError', ...
                'Error finding triangulation (barycentric) coordinates.' ...
            );
        end
        tsm = xff('new:tsm');
        tsm.NrOfTargetVertices = size(icov, 1);
        tsm.NrOfSourceVertices = size(srfv, 1);
        tsm.NrOfSourceTriangles = size(srft, 1);
        tsm.SourceTriangleOfTarget = t;
        tsm.TriangleEdgeLengths = vl;
        tsm.SaveAs(sprintf('%s_ICO%d.tsm', srff(1:end-4), size(icov, 1)));
        tsms{sc} = tsm.FilenameOnDisk;
        sph = srf.ApplyTSM(tsm, ico);
        sph.SaveAs(sprintf('%s_ICO%d.srf', srff(1:end-4), size(icov, 1)));
        srfs{sc} = sph.FilenameOnDisk;
        spd = smp.ApplyTSM(tsm, srf);
        spm.Map = smp.Map;
        for mc = 1:numel(spm.Map)
            spm.Map(mc).SMPData = spd(:, mc);
        end
        spm.NameOfOriginalSRF = strrep(srff, '_INFL_SPHERE_DC', '');
        spm.SaveAs(sprintf('%s_ICO%d_CURVATURE.smp', ...
            strrep(srff(1:end-4), '_INFL_SPHERE_DC', ''), size(icov, 1)));
        smps{sc} = spm.FilenameOnDisk;
        sph.ClearObject;
        smp.ClearObject;
        srf.ClearObject;
        tsm.ClearObject;
    end
    spm.ClearObject;
    ico.ClearObject;
end

% create BVX file for coordinate, curvature, trianglevertex data storage
hxr = hxdouble(randn(1, 2));
bvx = sprintf('%s/bvx_%s.bvx', tempdir, hxr([5:14, 21:30]));
bvxok = bvxcreatefile(bvx);
bvxo2 = bvxaddvartofile(bvx, ...
    {'nsrf', 'alstate', 'icomesh', 'nrvert', ...
     'temp_c', 'temp_v', 'temp_t', 'temp_m'}, ...
    {nsrf, zeros(8, max(80, nsrf)), double(opts.icomesh), opts.nrvert, ...
     ts_c, ts_v, ts_t, ts_m});
if ~bvxok || ...
   ~bvxo2
    error( ...
        'neuroelf:BVXError', ...
        'Error creating temporary storage file.' ...
    );
end

% loop over given surfaces
for sc = 1:nsrf
    try
        srf = [];
        smp = [];
        srf = xff(srfs{sc});
        smp = xff(smps{sc});
        if numel(srf) ~= 1 || ...
           ~isxff(srf, 'srf') || ...
            numel(smp) ~= 1 || ...
           ~isxff(smp, 'smp') || ...
            srf.NrOfVertices ~= smp.NrOfVertices || ...
            smp.NrOfMaps ~= 4
            error('BADOBJ');
        end
    catch ne_eo;
        if ~isempty(srf)
            srf.ClearObject;
        end
        if ~isempty(smp)
            smp.ClearObject;
        end
        delete(bvx);
        error( ...
            'neuroelf:xffError', ...
            'Error loading SRF/SMP %d (%s).', ...
            sc, ne_eo.message ...
        );
    end
    curv = [smp.Map(1).SMPData(:), smp.Map(2).SMPData(:), ...
            smp.Map(3).SMPData(:), smp.Map(4).SMPData(:)];
    bvxok = bvxaddvartofile(bvx, ...
        {sprintf('s%04d_c', sc), sprintf('s%04d_v', sc), ...
         sprintf('s%04d_t', sc), sprintf('s%04d_m', sc)}, ...
        {srf.MeshCenter, srf.VertexCoordinate, srf.TriangleVertex, curv});
    srf.ClearObject;
    smp.ClearObject;
    if ~bvxok
        delete(bvx);
        error( ...
            'neuroelf:BVXError', ...
            'Error adding information for SRF %d.', ...
            sc ...
        );
    end
end

% now call the actual alignment based on the BVX file
srfbvxalign(bvx);

% finally we can delete the bvx, ~whew~
delete(bvx);
