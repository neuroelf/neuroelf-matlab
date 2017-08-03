function [tc, tcf, tcv, tr] = mdm_VOITimeCourses(xo, voi, opts)
% MDM::VOITimeCourses  - extract time courses for referenced VTCs
%
% FORMAT:       [tc, tcf, tcv, tr] = mdm.VOITimeCourses(voi [, options])
%
% Input fields:
%
%       voi         VOI object
%       options     optional 1x1 struct with fields
%        .pbar      1x1 progress bar (xfigure::progress or xprogress)
%        .subsel    cell array with subject IDs to work on
%        .subvois   subject specific VOIs, either 'sub_', '_sub', {'voi'}
%        .trans     either of {'none'}, 'psc', 'z'
%        .unique    flag, only extract unique functional voxels (false)
%        .voisel    cell array with sub-VOI selection to use
%
% Output fields:
%
%       tc          Sx1 cell array with TxV time courses (S = NrOfStudies)
%       tcf         Sx1 cell array with time-course file names
%       tcv         Sx1 cell array with 1xV VOI indices used
%
% Using: findfirst, multimatch, psctrans, ztrans.

% Version:  v1.1
% Build:    16012715
% Date:     Jan-27 2016, 3:52 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

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

% neuroelf library and factory
global ne_methods xffsngl;
findfirst = ne_methods.findfirst;
psctrans  = ne_methods.psctrans;
ztrans    = ne_methods.ztrans;

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'mdm') || ...
    numel(voi) ~= 1 || ~xffisobject(voi, true, {'poi', 'voi'})
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
voitype = lower(voi.S.Extensions{1}(1));
cfs = struct('autofind', true, 'silent', true);
try
    mdm_CheckFiles(xo, cfs);
catch xfferror
    rethrow(xfferror);
end
bc = xo.C;
if (voitype == 'v' && (~strcmpi(bc.TypeOfFunctionalData, 'vtc') || size(bc.XTC_RTC, 2) ~= 2)) || ...
   (voitype == 'p' && (~strcmpi(bc.TypeOfFunctionalData, 'mtc')))
    error('neuroelf:xff:badArgument', ...
        'MDM must be VTC-based for VOI calls and MTC-based for POI calls.');
end

% check options
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'pbar') || numel(opts.pbar) ~= 1 || ...
   (~isxfigure(opts.pbar) && ~isa(opts.pbar, 'xprogress'))
    pbar = [];
else
    pbar = opts.pbar;
    pbarvis = pbar.Visible;
end
if ~isfield(opts, 'subsel') || ~iscell(opts.subsel) || isempty(opts.subsel)
    opts.subsel = mdm_Subjects(xo);
else
    opts.subsel = opts.subsel(:);
    try
        ssm = ne_methods.multimatch(opts.subsel, mdm_Subjects(xo));
    catch xfferror
        rethrow(xfferror);
    end
    if any(ssm < 1)
        error('neuroelf:xff:badArgument', 'Invalid subject ID in selection.');
    end
end
subids = mdm_Subjects(xo, true);
subseli = ne_methods.multimatch(subids, opts.subsel);
nsubsel = sum(subseli > 0);
if nsubsel == 0
    error('neuroelf:xff:badArgument', 'No subjects selected.');
end
if ~isfield(opts, 'subvois') || ~ischar(opts.subvois) || isempty(opts.subvois) || ...
   ~any(lower(opts.subvois(1)) == '_s')
    opts.subvois = 'v';
else
    opts.subvois = lower(opts.subvois(1));
end
if ~isfield(opts, 'trans') || ~ischar(opts.trans) || isempty(opts.trans) || ...
   ~any('npz' == lower(opts.trans(1)))
    opts.trans = 'n';
    if bc.PSCTransformation > 0
        opts.trans = 'p';
    end
    if bc.zTransformation > 0
        opts.trans = 'z';
    end
else
    opts.trans = lower(opts.trans(1));
end
if ~isfield(opts, 'unique') || ~islogical(opts.unique) || numel(opts.unique) ~= 1
    opts.unique = false;
end
if voitype == 'v'
    vois = voi_VOINames(voi);
else
    vois = poi_POINames(voi);
end
if ~isfield(opts, 'voisel') || ~iscell(opts.voisel) || isempty(opts.voisel)
    switch (opts.subvois)
        case {'_'}
            svois = unique(lower(regexprep(vois, '_[^_]+$', '')));
            for vc = numel(svois):-1:1
                vf = true;
                for sc = 1:numel(opts.subsel)
                    if ~any(strcmpi(vois, [svois{vc} '_' opts.subsel{sc}]))
                        vf = false;
                        break;
                    end
                end
                if ~vf
                    svois(vc) = [];
                    continue;
                end
            end
            if isempty(svois)
                error('neuroelf:xff:badArgument', 'No VOI is available for all subjects.');
            end
            opts.voisel = svois;
        case {'s'}
            svois = unique(lower(regexprep(vois, '^[^_]+_', '')));
            for vc = numel(svois):-1:1
                vf = true;
                for sc = 1:numel(opts.subsel)
                    if ~any(strcmpi(vois, [opts.subsel{sc} '_' svois{vc}]))
                        vf = false;
                        break;
                    end
                end
                if ~vf
                    svois(vc) = [];
                    continue;
                end
            end
            if isempty(svois)
                error('neuroelf:xff:badArgument', 'No VOI is available for all subjects.');
            end
            opts.voisel = svois;
        case {'v'}
            opts.voisel = vois;
    end
else
    opts.voisel = opts.voisel(:);
    try
        switch (opts.subvois)
            case {'_'}
                for vc = 1:numel(opts.voisel)
                    for sc = 1:numel(opts.subsel)
                        if ~any(strcmpi(vois, [opts.voisel{vc} '_' opts.subsel{sc}]))
                            error('neuroelf:xff:badArgument', ...
                                'VOI %s not available for subject %s.', ...
                                opts.voisel{vc}, opts.subsel{sc});
                        end
                    end
                end
            case {'s'}
                for vc = 1:numel(opts.voisel)
                    for sc = 1:numel(opts.subsel)
                        if ~any(strcmpi(vois, [opts.subsel{sc} '_' opts.voisel{vc}]))
                            error('neuroelf:xff:badArgument', ...
                                'VOI %s not available for subject %s.', ...
                                opts.voisel{vc}, opts.subsel{sc});
                        end
                    end
                end
            case {'v'}
                for vc = 1:numel(opts.voisel)
                    if ~any(strcmpi(vois, opts.voisel{vc}))
                        error('neuroelf:xff:badArgument', ...
                            'VOI %s not available.', opts.voisel{vc});
                    end
                end
        end
    catch xfferror
        rethrow(xfferror);
    end
end

% generate output arguments
tc = cell(nsubsel, 1);
if nargout > 1
    tcf = cell(nsubsel, 1);
    if nargout > 2
        tcv = cell(nsubsel, 1);
        if nargout > 3
            tr = zeros(nsubsel, 1);
        end
    end
end

% set transio for VTCs
vtcidx = xffsngl.EXT.vtc{2};
tios = xffsngl.BFF(vtcidx).TransIOSize;
xffsngl.BFF(vtcidx).TransIOSize = 8192;
xffsngl.FF.bff(vtcidx).TransIOSize = 8192;

% copy VOI object and get content
cvoi = aft_CopyObject(voi);
cvoic = cvoi.C;
ovoic = cvoic;
if voitype == 'v'
    cvoin = {cvoic.VOI.Name};
else
    cvoin = {cvoic.POI.Name};
end

% VOI indices
vid = zeros(1, numel(opts.voisel));
cvoic.NrOfVOIs = numel(vid);

% VOI pattern (without subject ID in name)
if opts.subvois == 'v'
    for vc = 1:numel(opts.voisel)
        vid(vc) = findfirst(strcmpi(cvoin, opts.voisel{vc}));
    end
    if voitype == 'v'
        cvoic.VOI = ovoic.VOI(vid);
    else
        cvoic.POI = ovoic.POI(vid);
    end
    cvoi.C = cvoic;
end

% to-clear objects
tco = {cvoi, []};

% for progress
if ~isempty(pbar)

    % initialize progress
    if voitype == 'v'
        pbar.Progress(0, 'Extracting VOI time courses...');
    else
        pbar.Progress(0, 'Extracting POI time courses...');
    end
    pbar.Visible = 'on';
    drawnow;
end

% big TRY/CATCH
try

    % iterate over subsel
    tcfi = 1;
    for sc = 1:numel(opts.subsel)

        % get required time course files
        tcfs = bc.XTC_RTC(subseli == sc, 1);

        % iterate over time courses
        for fc = 1:numel(tcfs)

            % store name
            if nargout > 1
                tcf{tcfi} = tcfs{fc};
            end

            % open VTC
            tco{2} = xff(tcfs{fc});

            % get TR if necessary
            if nargout > 3
                tcoc = tco{2}.C;
                tr(tcfi) = tcoc.TR;
            end

            % depending on access
            switch (opts.subvois)

                % VOI_SUB pattern
                case {'_'}

                    % iterate over selected VOIs
                    for vc = 1:numel(opts.voisel)
                        vid(vc) = findfirst(strcmpi(cvoin, ...
                            [opts.subsel{sc} '_' opts.voisel{vc}]));
                    end

                % SUB_VOI pattern
                case {'s'}
                    for vc = 1:numel(opts.voisel)
                        vid(vc) = findfirst(strcmpi(cvoin, ...
                            [opts.voisel{vc} '_' opts.subsel{sc}]));
                    end
            end

            % patch copied VOI
            if opts.subvois ~= 'v'
                if voitype == 'v'
                    cvoic.VOI = ovoic.VOI(vid);
                else
                    cvoic.POI = ovoic.POI(vid);
                end
                cvoi.C = cvoic;
            end

            % extract time courses
            tc{tcfi} = aft_VOITimeCourse(tco{2}, cvoi, ...
                struct('weights', 1 - double(opts.unique)));
            aft_ClearObject(tco{2});
            tco{2} = [];

            % transformation
            switch (opts.trans)
                case {'p'}
                    tc{tcfi} = psctrans(tc{tcfi}, 1);
                case {'z'}
                    tc{tcfi} = ztrans(tc{tcfi}, 1);
            end

            % VOI names
            if nargout > 2
                tcv{tcfi} = vid;
            end

            % keep UI alive
            if ~isempty(pbar)
                pbar.Progress(tcfi / nsubsel);
            end
            drawnow;

            % increase counter
            tcfi = tcfi + 1;
        end
    end
catch xfferror

    % re-set transiosize
    xffsngl.BFF(vtcidx).TransIOSize = tios;
    xffsngl.FF.bff(vtcidx).TransIOSize = tios;

    % free copied VOI/VTC object
    clearxffobjects(tco);

    % bail out
    rethrow(xfferror);
end

% hide progress bar?
if ~isempty(pbar)
    pbar.Visible = pbarvis;
    drawnow;
end

% re-set transiosize
xffsngl.BFF(vtcidx).TransIOSize = tios;
xffsngl.FF.bff(vtcidx).TransIOSize = tios;

% free copied VOI object as well as VTC object
clearxffobjects(tco);
