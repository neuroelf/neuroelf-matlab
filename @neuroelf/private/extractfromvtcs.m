function [xtc, xstd, opts] = extractfromvtcs(files, xfiles, opts)
% extractfromvtcs  - extract time courses from VTCs
%
% FORMAT:       [xtc, xstd, opts] = extractfromvtcs(files, xfiles [, opts])
%
% Input fields:
%
%       files       either SxR cell array or folder name + pattern
%       xfiles      filenames (or pattern) of either MSKs, SRFs, or VOIs
%       opts        optional settings
%        .savetcs   filename pattern to save timecourses ('%S_%R_%V.sdm')
%        .stdtrans  either 'none', {'psc'}, or 'z'
%        .tctrans   either {'none'}, 'psc', or 'z'
%
% Output fields:
%
%       xtc         SxR cell array (subjects-by-runs) with TxV time courses
%       xstd        SxRxV standard deviations
%       opts        used options with additional fields
%        .runids    Rx1 cell array of run IDs
%        .subids    Sx1 cell array of subject IDs
%        .voiids    Vx1 cell array of VOI IDs

% Version:  v0.9c
% Build:    12020811
% Date:     Nov-29 2011, 12:25 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, Jochen Weber
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
    isempty(files) || ...
    isempty(xfiles) || ...
    ndims(files) > 2 || ...
    ndims(xfiles) > 2 || ...
   (~iscell(files) && ...
    ~ischar(files)) || ...
   (ischar(files) && ...
    ~any(files(:) == '*')) || ...
   (~iscell(xfiles) && ...
    ~ischar(xfiles))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin < 3 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'savetcs') || ...
   ~ischar(opts.savetcs)
    opts.savetcs = '%S_%R_%V.sdm';
else
    opts.savetcs = opts.savetcs(:)';
    if ~isempty(opts.savetcs) && ...
       (numel(opts.savetcs) < 10 || ...
        isempty(opts.savetcs, '%S') || ...
        isempty(opts.savetcs, '%R') || ...
        isempty(opts.savetcs, '%V') || ...
        ~any(strcmpi(opts.savetcs(end-3:end), {'.mat', '.sdm', '.txt'})))
        error( ...
            'neuroelf:BadArgument', ...
            'Bad or missing savetcs option.' ...
        );
    end
end
if ~isempty(opts.savetcs)
    stext = lower(opts.savetcs(end-2));
else
    stext = '';
end
if ~isfield(opts, 'stdtrans') || ...
   ~ischar(opts.stdtrans) || ...
    isempty(opts.stdtrans) || ...
   ~any(lower(opts.stdtrans(1)) == 'npz')
    opts.stdtrans = 'p';
else
    opts.stdtrans = lower(opts.stdtrans(1));
end
if ~isfield(opts, 'tctrans') || ...
   ~ischar(opts.tctrans) || ...
    isempty(opts.tctrans) || ...
   ~any(lower(opts.tctrans(1)) == 'npz')
    opts.tctrans = 'n';
else
    opts.tctrans = lower(opts.tctrans(1));
end

% parse files
if ischar(files)

    % for a pattern
    try

        % look for files
        filesp = fileparts(files);
        if any(filesp == '*')
            files = findfiles(files);
        else
            files = findfiles(files, 'depth=1');
        end

    % deal with errors
    catch ne_eo;
        rethrow(ne_eo);
    end

% for cell arrays
else

    % vectorize
    files = files(:);
end

% then parse into subjects and runs
subids = files;
runids = files;
for fc = 1:numel(files)
    [fdir, fname] = fileparts(files{fc});
    fparts = splittocellc(fname, '_');
    if numel(fparts) < 2
        error( ...
            'neuroelf:BadArgument', ...
            'Invalid VTC filename: %s.', ...
            files{fc} ...
        );
    end
    subids{fc} = fparts{1};
    runids{fc} = fparts{2};
end

% then create SxR array
opts.runids = uunion(runids(end:-1:1), runids);
opts.runids = opts.runids(end:-1:1);
opts.subids = uunion(subids(end:-1:1), subids);
opts.subids = opts.subids(end:-1:1);
farray = cell(numel(opts.subids), numel(opts.runids));

% for each unique subject ID
for sc = 1:numel(opts.subids)

    % get entries
    subruns = find(strcmp(opts.subids{sc}, subids));

    % and divide into runs
    subruni = multimatch(runids(subruns), opts.runids);
    if numel(subruni) ~= numel(unique(subruni))
        error( ...
            'neuroelf:BadArgument', ...
            'Invalid run IDs for VTCs of subject %s.', ...
            opts.subids{sc} ...
        );
    end

    % then store filenames
    farray(sc, subruni) = files(subruns)';
end

% now we need to parse the extraction files
if ischar(xfiles)

    % for a pattern
    try

        % look for files
        filesp = fileparts(xfiles);
        if any(filesp == '*')
            xfiles = findfiles(xfiles);
        else
            xfiles = findfiles(xfiles, 'depth=1');
        end

    % deal with errors
    catch ne_eo;
        rethrow(ne_eo);
    end

% for cell arrays
else

    % vectorize
    xfiles = xfiles(:);
end

% types of files
msk = false;
msklay = [];
srf = false;
voi = false;
xsubids = xfiles;
xvoiids = xfiles;
xobjs = xfiles;
for fc = 1:numel(xfiles)

    % inspect filename
    [fdir, fname, fext] = fileparts(xfiles{fc});
    if isempty(fext)
        error( ...
            'neuroelf:BadArgument', ...
            'Extract files must have an extension.' ...
        );
    end
    fext = lower(fext(2:end));
    if ~any(strcmp(fext, {'msk'; 'srf'; 'voi'}))
        error( ...
            'neuroelf:BadArgument', ...
            'Extract files must have either a .msk, .srf, or .voi extension.' ...
        );
    end
    try
        switch (fext)
            case {'msk'}
                msk = true;
            case {'srf'}
                srf = true;
            case {'voi'}
                voi = true;
                xobjs{fc} = xff(xfiles{fc});
                if numel(xobjs{fc}) ~= 1 || ...
                   ~isxff(xobjs{fc}, 'voi')
                    error( ...
                        'neuroelf:BadArgument', ...
                        'Not a VOI file: %s.', ...
                        xfiles{fc} ...
                    );
                end
        end
    catch ne_eo;
        clearxffobjects(xobjs);
        rethrow(ne_eo);
    end

    % check for mixing
    if (double(msk) + double(srf) + double(voi)) ~= 1
        error( ...
            'neuroelf:BadArgument', ...
            'This function does not support mixing extraction filetypes.' ...
        );
    end

    % get subject ID
    fparts = splittocellc(fname, '_');

    % test against existing list
    if ~any(strcmp(opts.subids, fparts{1}))

        % no match and not VOI
        if ~voi

            % error out
            clearxffobjects(xobjs);
            error( ...
                'neuroelf:BadArgument', ...
                'Extract files must match to VTC subject IDs.' ...
            );
        end
    end

    % for VOIs
    if voi

        % make sure VOI names include subject ID
        voinames = xobjs{fc}.VOINames;
        if ~any(strcmp(opts.subids, fparts{1}))
            for vc = numel(voinames):-1:1
                if isempty(xobjs{fc}.VOI(vc).Voxels)
                    xobjs{fc}.VOI(vc) = [];
                    continue;
                end
                vparts = splittocellc(voinames{vc}, '_');
                if ~any(strcmp(opts.subids, vparts{1}))
                    clearxffobjects(xobjs);
                    error( ...
                        'neuroelf:BadArgument', ...
                        'VOI names of files without subject ID must contain a subject ID.' ...
                    );
                end
            end
        else
            for vc = numel(voinames):-1:1
                if isempty(xobjs{fc}.VOI(vc).Voxels)
                    xobjs{fc}.VOI(vc) = [];
                    continue;
                end
                vparts = splittocellc(voinames{vc}, '_');
                if ~strcmp(vparts{1}, fparts{1})
                    if any(strcmp(opts.subids, vparts{1}))
                        error( ...
                            'neuroelf:BadArgument', ...
                            'Genuine VOI names must not match another subject''s ID.' ...
                        );
                    else
                        xobjs{fc}.VOI(vc).Name = sprintf('%s_%s', fparts{1}, voinames{vc});
                    end
                end
            end
        end

    % otherwise
    else

        % require at least two parts
        if numel(fparts) < 2
            clearxffobjects(xobjs);
            error( ...
                'neuroelf:BadArgument', ...
                'MSK/SRF filenames must contain a VOI ID as second particle.' ...
            );
        end

        % make sure file exists and is readable
        try
            xobjs{fc} = xff(xfiles{fc});
            if numel(xobjs{fc}) ~= 1 || ...
               ~isxff(xobjs{fc}, fext)
                error( ...
                    'neuroelf:BadArgument', ...
                    'Invalid extraction file: %s.', ...
                    xfiles{fc} ...
                );
            end

            % for masks
            if msk

                % test layout
                if isempty(msklay)
                    msklay = xobjs{fc}.Layout;
                else
                    if ~isequal(msklay, xobjs{fc}.Layout)
                        error( ...
                            'neuroelf:BadArgument', ...
                            'MSK layouts must match.' ...
                        );
                    end
                end
            end

            % but don't accumulate memory!
            xobjs{fc}.ClearObject;
            xobjs{fc} = xfiles{fc};
        catch ne_eo;
            clearxffobjects(xobjs);
            rethrow(ne_eo);
        end
    end

    % depending on type
    if ~voi

        % store subject ID
        xsubids{fc} = fparts{1};
        xvoiids{fc} = fparts{2};
    else

        % get list
        voinames = xobjs{fc}.VOINames;
        svoiids = voinames;
        for vc = 1:numel(voinames)
            vparts = splittocellc(voinames{vc}, '_');
            svoiids{vc} = vparts{2};
        end
        xvoiids{fc} = svoiids;
    end
end

% for VOIs
if voi

    % join VOI ID list
    jvoiids = cat(1, xvoiids{:});
    opts.voiids = uunion(jvoiids(end:-1:1), jvoiids);

% for other types
else
    opts.voiids = uunion(xvoiids(end:-1:1), xvoiids);
end
opts.voiids = opts.voiids(end:-1:1);
nvois = numel(opts.voiids);

% now we can already create the outputs
xtc = cell(size(farray));
xstd = NaN .* ones(numel(farray), nvois);

% then try to load each (using transio)
objs = farray;
try
    for fc = 1:numel(objs)

        % if this subject/run combination isn't available
        if isempty(objs{fc})
            continue;
        end

        % otherwise get object handle
        objs{fc} = xff(objs{fc}, 't');

        % test object
        if numel(objs{fc}) ~= 1 || ...
           ~isxff(objs{fc}, 'vtc')
            error( ...
                'neuroelf:BadArgument', ...
                'Not a VTC file: %s.', ...
                files{fc} ...
            );
        end

        % for msk access, test layout
        if msk
            vtclay = objs{fc}.Layout;
            if any(msklay(5:13) ~= vtclay(5:13))
                error( ...
                    'neuroelf:BadArgument', ...
                    'MSK and VTC layout must match.' ...
                );
            end
        end

        % and create time course array within
        xtc{fc} = nan .* ones(objs{fc}.NrOfVolumes, nvois);
    end
catch ne_eo;
    clearxffobjects(objs);
    clearxffobjects(xobjs);
    rethrow(ne_eo);
end

% extract from time courses
for fc = 1:numel(xobjs)

    % for VOIs
    if voi

        % get VOI names
        voinames = xobjs{fc}.VOINames;

        % for each VOI
        for vc = 1:numel(voinames)

            % get subject and VOI ID
            vparts = splittocellc(voinames{vc}, '_');

            % match to list
            sid = find(strcmp(vparts{1}, opts.subids));
            vid = find(strcmp(vparts{2}, opts.voiids));

            % iterate over runs for that subject
            for rc = 1:size(farray, 2)
                if ~isempty(objs{sid, rc})
                    disp(sprintf('Extracting (%s, %s, %s)', ...
                        opts.subids{sid}, opts.runids{rc}, opts.voiids{vid}));
                    pause(0.001);
                    xtc{sid, rc}(:, vid) = ...
                        objs{sid, rc}.VOITimeCourse(xobjs{fc}, struct('voisel', vc));
                end
            end
        end

    % for other objects
    else

        % get subject and VOI ID
        sid = find(strcmp(xsubids{fc}, opts.subids));
        vid = find(strcmp(xvoiids{fc}, opts.voiids));

        % load object
        xobjs{fc} = xff(xobjs{fc});

        % iterate over runs for that subject
        for rc = 1:size(farray, 2)

            % extract time course
            if ~isempty(objs{sid, rc})

                % info
                disp(sprintf('Extracting (%s, %s, %s)', ...
                    opts.subids{sid}, opts.runids{rc}, opts.voiids{vid}));
                pause(0.001);

                % for MSK
                if msk
                    xtc{sid, rc}(:, vid) = meannoinfnan( ...
                        objs{sid, rc}.VTCData(:, find(xobjs{fc}.Mask(:))), 2, true);

                % for SRF
                else
                    xsdm = objs{sid, rc}.MeanSRFTimeCourse(xobjs{fc}, ...
                        struct('remmean', false));
                    xtc{sid, rc}(:, vid) = xsdm.SDMMatrix;
                    xsdm.ClearObject;
                end
            end
        end

        % clear object
        xobjs{fc}.ClearObject;
        xobjs{fc} = xfiles{fc};
    end
end

% clear objects
clearxffobjects(xobjs);
clearxffobjects(objs);

% computations
for fc = 1:numel(xtc)

    % nothing for empty arrays
    if isempty(xtc{fc})
        continue;
    end

    % what type of transformation for STD computation
    switch (opts.stdtrans)
        case {'n'}
            xstd(fc, :) = std(xtc{fc}, 0, 1);
        case {'p'}
            xstd(fc, :) = std(psctrans(xtc{fc}), 0, 1);
        case {'z'}
            xstd(fc, :) = std(ztrans(xtc{fc}), 0, 1);
    end

    % what type of transformation
    switch (opts.tctrans)
        case {'n'}
            xtc{fc} = xtc{fc};
        case {'p'}
            xtc{fc} = psctrans(xtc{fc});
        case {'z'}
            xtc{fc} = ztrans(xtc{fc});
    end
end

% reshape xstd
xstd = reshape(xstd, [size(farray), nvois]);

% save outputs
switch (stext)

    % MAT or text files
    case {'m', 't'}

        % for text files
        if stext == 't'
            ascflag = {'-ascii'};

        % but not for MAT files
        else
            ascflag = {};
        end

        % iterate over data
        for sc = 1:size(xtc, 1)
            for rc = 1:size(xtc, 2)
                if isempty(xtc{sc, rc})
                    continue;
                end
                for vc = 1:nvois
                    if isnan(xstd(fc, vc))
                        continue;
                    end

                    % generate name
                    fname = strrep(strrep(strrep(opts.savetcs, ...
                        '%S', opts.subids{sc}), '%R', opts.runids{rc}), ...
                        '%V', opts.voiids{vc});

                    % and save time course
                    tc = xtc{sc, rc}(:, vc);
                    save(fname, 'tc', ascflag{:});
                end
            end
        end

    % SDM files
    case {'s'}

        % generate SDM
        sdm = xff('new:sdm');
        sdm.NrOfPredictors = 1;
        sdm.IncludesConstant = 0;
        sdm.FirstConfoundPredictor = 2;
        sdm.PredictorColors = [128, 128, 128];

        % iterate over data
        for sc = 1:size(xtc, 1)
            for rc = 1:size(xtc, 2)
                if isempty(xtc{sc, rc})
                    continue;
                end

                % adapt SDM (part 1)
                sdm.NrOfDataPoints = size(xtc{sc, rc}, 1);

                % iterate over VOIs
                for vc = 1:nvois
                    if isnan(xstd(sc, rc, vc))
                        continue;
                    end

                    % generate name
                    fname = strrep(strrep(strrep(opts.savetcs, ...
                        '%S', opts.subids{sc}), '%R', opts.runids{rc}), ...
                        '%V', opts.voiids{vc});

                    % adapt SDM
                    sdm.PredictorNames{1} = sprintf('%s_%s_%s', ...
                        opts.subids{sc}, opts.runids{rc}, opts.voiids{vc});
                    sdm.SDMMatrix = xtc{sc, rc}(:, vc);
                    sdm.SaveAs(fname);
                end
            end
        end

        % clear object
        sdm.ClearObject;
end
