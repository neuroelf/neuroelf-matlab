function sdms = testjittershifts(mdm, remodisis, opts)
% testjittershifts  - test different shifts in jitters (remodeled ISIs)
%
% FORMAT:       sdms = testjittershifts(mdm, remodisis [, opts))
%
% Input fields:
%
%       mdm         MDM object (or filename)
%       remodisis   remodisis flag for MDM::SDMs call
%       opts        optional settings
%        .cons      contrasts (cell array) for which to compute cbcov
%        .ffx       within-subjects, across-runs FFX models (default: false)
%        .prtr      if set, create VTCs with altered TR (will be deleted)
%        .sdmcut    cut SDMs to portion with data in conditions (true)
%
% Output fields:
%
%       sdms        SDM structures with added fields for VIF and cBcov

% Version:  v0.9d
% Build:    14061709
% Date:     Jun-17 2014, 9:47 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2013, 2014, Jochen Weber
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
   (~isxff(mdm, 'mdm') && ...
    (~ischar(mdm) || ...
     isempty(mdm) || ...
     exist(mdm(:)', 'file') ~= 2 || ...
     isempty(regexpi(mdm(:)', '\.mdm$')))) || ...
   ~iscell(remodisis) || ...
    ndims(remodisis) ~= 2 || ...
   ~any(size(remodisis, 2) == [2, 3])
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
if ~isfield(opts, 'cons') || ...
   ~iscell(opts.cons) || ...
    isempty(opts.cons)
    opts.cons = {};
else
    opts.cons = opts.cons(:);
    for cc = numel(opts.cons):-1:1
        if ~ischar(opts.cons{cc}) || ...
            isempty(opts.cons{cc})
            opts.cons(cc) = [];
        end
    end
end
if ~isfield(opts, 'ffx') || ...
   ~islogical(opts.ffx) || ...
    numel(opts.ffx) ~= 1
    opts.ffx = false;
end
prtr = [];
if isfield(opts, 'prtr') && ...
    isa(opts.prtr, 'double') && ...
    numel(opts.prtr) == 1 && ...
   ~isinf(opts.prtr) && ...
   ~isnan(opts.prtr) && ...
    opts.prtr >= 50
    prtr = opts.prtr;
end
if ~isfield(opts, 'sdmcut') || ...
   ~islogical(opts.sdmcut) || ...
    numel(opts.sdmcut) ~= 1
    opts.sdmcut = true;
end

% load MDM if needed
mdmloaded = false;
if ischar(mdm)
    try
        mdm = xff(mdm(:)');
    catch ne_eo;
        rethrow(ne_eo);
    end
    mdmloaded = true;
end

% all design files need to be PRT
if any(cellfun('isempty', regexpi(mdm.XTC_RTC(:, 2), '\.prt$')))
    if mdmloaded
        mdm.ClearObject;
    end
    error( ...
        'neuroelf:BadArgument', ...
        'MDM must be PRT-based.' ...
    );
end
studies = mdm.XTC_RTC;
numstudies = size(studies, 1);

% for fixed effects
if opts.ffx

    % get subject IDs (short and full)
    subjects = mdm.Subjects;
    runsubids = mdm.Subjects(true);
end

% new TR requires to recreate VTCs
if ~isempty(prtr)

    % generate VTC with 2x2x2 voxels
    vtc = xff('new:vtc');
    vtc.TR = prtr;
    vtc.Resolution = 3;
    vtc.XStart = 57;
    vtc.XEnd = 63;
    vtc.YStart = 52;
    vtc.YEnd = 58;
    vtc.ZStart = 59;
    vtc.ZEnd = 65;
    vtc.DataType = 2;
    vtc.FileVersion = 3;
    vtcfiles = cell(numstudies, 1);

    % iterate over studies
    for sc = 1:numstudies
        try
            prt = xff(studies{sc, 2});
            if ~strcmpi(prt.ResolutionOfTime, 'msec')
                prt.ClearObject;
                error( ...
                    'neuroelf:BadObject', ...
                    'PRTs must be millisecond-based.' ...
                );
            end
            if ~isempty(remodisis)
                prt.RemodelISIs(remodisis);
                prt.CleanUp;
            end
            prtc = prt.Cond;
            prt.ClearObject;
            prto = cat(1, prtc.OnOffsets);
            lastonset = max(prto(:)) + 15000;
            numvols = 1 + ceil(lastonset / prtr);

            % generate new VTCData
            vtc.VTCData = single(100 + 5 .* randn([numvols, 2, 2, 2]));

            % save as
            vtcfile = [tempname '.vtc'];
            vtc.SaveAs(vtcfile);
            vtcfiles{sc} = vtcfile;
        catch ne_eo;
            ne_eoo = ne_eo;
            if mdmloaded
                mdm.ClearObject;
            end
            for ssc = 1:sc
                try
                    delete(vtcfiles{ssc});
                catch ne_eo;
                    warning(ne_eo.message);
                end
            end
            rethrow(ne_eoo);
        end
    end

    % update file list
    studies(:, 1) = vtcfiles;
    mdm.XTC_RTC = studies;
end

% augment opts
opts.asstruct = true;
if ~isempty(remodisis)
    opts.remodisis = remodisis;
end

% get SDMs
sdms = mdm.SDMs(opts);

% delete temporary VTCs
if ~isempty(prtr)
    for sc = 1:numel(vtcfiles)
        try
            delete(vtcfiles{sc});
        catch ne_eo;
            warning(ne_eo.message);
        end
    end
end

% unload MDM
if mdmloaded
    mdm.ClearObject;
end

% truncate SDMs
if opts.sdmcut
    for sc = 1:numel(sdms)
        firstallzero = findfirst(any(sdms{sc}.SDMMatrix(:, 1:sdms{sc}.FirstConfoundPredictor-1) ~= 0, 2), -1);
        if firstallzero < size(sdms{sc}.SDMMatrix, 1)
            sdms{sc}.SDMMatrix(firstallzero+1:end, :) = [];
        end
    end
end

% for FFX
if opts.ffx

    % generate ffx SDMs array
    ffxsdms = cell(numel(subjects), 1);

    % iterate over subjects
    for sc = 1:numel(subjects)

        % find matching runs
        runmatch = find(strcmpi(runsubids, subjects{sc}));

        % iterate over found runs
        for rc = 1:numel(runmatch)

            % for first run
            rmc = runmatch(rc);
            if rc == 1

                % simply copy SDM
                ffxsdms{sc} = sdms{rmc};
                X = ffxsdms{sc}.SDMMatrix;
                Xp = ffxsdms{sc}.PredictorNames(:);
                Xc = ffxsdms{sc}.PredictorColors;
                Xf = ffxsdms{sc}.FirstConfoundPredictor;
                continue;
            end

            % get size of second design
            sX = sdms{rmc}.SDMMatrix;
            fsz = size(X);
            dsz = size(sX);
            fcp = sdms{rmc}.FirstConfoundPredictor;

            % then match then regressors (of interest) by name
            rnames1 = Xp(1:Xf-1);
            rnames2 = sdms{rmc}.PredictorNames(1:fcp-1);
            rnames2 = rnames2(:);
            rnamesm = multimatch(rnames2, rnames1);

            % unfound names
            if any(rnamesm == 0)

                % add to data, names and colors
                X = [X(:, 1:Xf-1), zeros(fsz(1), sum(rnamesm == 0)), X(:, Xf:end)];
                Xp = [rnames1; rnames2(rnamesm == 0); Xp(Xf:end)];
                Xc = [Xc(1:Xf-1, :); sdms{rmc}.PredictorColors(rnamesm == 0, :); Xc(Xf:end, :)];
                Xf = Xf + sum(rnamesm == 0);
                fsz = size(X);
                rnames1 = Xp(1:Xf-1);

                % rematch
                rnamesm = multimatch(rnames2, rnames1);
            end

            % extend in time dimension
            X(fsz(1)+1:fsz(1)+dsz(1), :) = 0;
            X(:, fsz(2)+1:fsz(2)+dsz(2)-(fcp-1)) = 0;

            % store data
            X(fsz(1)+1:end, rnamesm) = sX(:, 1:fcp-1);
            X(fsz(1)+1:end, fsz(2)+1:end) = sX(:, fcp:end);

            % store back
            if rc == numel(runmatch)
                ffxsdms{sc}.SDMMatrix = X;
                ffxsdms{sc}.PredictorNames = Xp(:)';
                ffxsdms{sc}.PredictorColors = Xc;
                ffxsdms{sc}.FirstConfoundPredictor = Xf;
            end
        end
    end

    % replace data
    sdms = ffxsdms;
end

% compute VIFs and inverse of design matrices for
for sc = 1:numel(sdms)
    sdms{sc}.VIF = vifactor(sdms{sc}.SDMMatrix);
    sdms{sc}.Bcov = inv(sdms{sc}.SDMMatrix' * sdms{sc}.SDMMatrix);
end

% compute contrast error factors
if ~isempty(opts.cons)
end
