function [sdm, sdms] = spmmat2sdm(spmmat, sdmfile, regoi)
% spmmat2sdm  - convert a SPM.mat file into BV's SDM file(s)
%
% FORMAT:       [sdm, sdms] = spmmat2sdm(spmmat, sdmfile [, regoi])
%
% Input fields:
%
%       spmmat      either SPM.mat filename or struct with fields
%        .SPM       1x1 struct containing loaded struct
%             - or -
%        .xX        design structure
%          .X       design matrix
%          .name    names of predictors
%       sdmfile     outfile filename
%       regoi       regressors of interest (all others put to confounds)
%
% Output fields:
%
%       sdm         SDM object of first session
%       sdms        cell array with SDM objects of all session
%
% See also spmmat2prt

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
if nargin < 2 || ...
   (~ischar(spmmat) && ...
    ~isstruct(spmmat)) || ...
    isempty(spmmat) || ...
   ~ischar(sdmfile) || ...
    isempty(sdmfile)
    error( ...
        'neuroelf:BadArgument', ...
        'Missing or bad argument. Try ''help %s''.', ...
        mfilename ...
    );
end
if nargin < 3 || ...
   ~iscell(regoi)
    regoi = {};
else
    regoi = regoi(:);
    for rc = numel(regoi):-1:1
        if ~ischar(regoi{rc}) || ...
            isempty(regoi{rc})
            regoi(rc) = [];
        else
            regoi{rc} = regoi{rc}(:)';
        end
    end
    regoi(strcmpi(regoi, 'constant')) = [];
    if ~isempty(regoi)
        regoi = uunion(regoi, {});
    end
end

% what is spmmat
if ischar(spmmat)
    try
        spmmat = load(spmmat(:)');
        spmmat.SPM;
    catch ne_eo;
        error( ...
            'neuroelf:BadArgument', ...
            'Invalid SPM.mat filename given (%s).', ...
            ne_eo.message ...
        );
    end
end
if numel(spmmat) ~= 1
    error( ...
        'neuroelf:BadArgument', ...
        'Bad dim/size of spmmat argument.' ...
    );
end
if isfield(spmmat, 'SPM')
    spmmat = spmmat.SPM;
end
if numel(spmmat) ~= 1 || ...
   ~isfield(spmmat, 'xX') || ...
   ~isstruct(spmmat.xX) || ...
    numel(spmmat.xX) ~= 1 || ...
   ~isfield(spmmat.xX, 'X') || ...
    isempty(spmmat.xX.X) || ...
   (~isa(spmmat.xX.X, 'double') && ...
    ~isa(spmmat.xX.X, 'single')) || ...
   ~isfield(spmmat.xX, 'name') || ...
   ~iscell(spmmat.xX.name) || ...
    numel(spmmat.xX.name) ~= size(spmmat.xX.X, 2)
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid SPM.mat structure given.' ...
    );
end
xX = spmmat.xX.X;
nampred = spmmat.xX.name(:);

% parse session numbers
numpred = numel(nampred);
sespred = zeros(1, numpred);
baspred = zeros(1, numpred);
for pc = 1:numpred
    if ~ischar(nampred{pc}) || ...
        numel(nampred{pc}) < 7 || ...
       ~strcmpi(nampred{pc}(1:3), 'sn(')
        error( ...
            'neuroelf:BadArgument', ...
            'Session not detectable from predictor names.' ...
        );
    end
    snnum = nampred{pc}(4:5);
    snnum(snnum == ')') = [];
    try
        snnum = str2double(snnum);
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        error( ...
            'neuroelf:BadArgument', ...
            'Session not detectable from predictor names.' ...
        );
    end
    sespred(pc) = snnum;
    if ~isempty(regexpi(nampred{pc}, '.*constant$'))
        baspred(pc) = 1;
    end
end
sessno = numel(unique(sespred));

% generate SDMs
[sdmabs{1:2}] = isabsolute(sdmfile(:)');
[sdmp{1:3}] = fileparts(sdmabs{2});
sdms = cell(1, sessno);
for sc = 1:sessno

    % get baseline column (to detect time indices)
    bsi = find((baspred > 0) & (sespred == sc));
    if numel(bsi) ~= 1
        clearxffobjects(sdms);
        error( ...
            'neuroelf:InternalError', ...
            'Error retrieving unique baseline predictor for session %d.', ...
            sc ...
        );
    end
    sti = find(xX(:, bsi) > 0);
    nti = numel(sti);

    % get predictor columns and names
    csi = find((baspred == 0) & (sespred == sc));
    csn = nampred(csi);
    for pc = 1:numel(csn)
        csn{pc} = regexprep(csn{pc}, ...
            '^sn\(\d+\)\s*([^\s\*]*)\*?.*$', '$1', 'preservecase');
    end

    % build SDM
    newsdm = xff('new:sdm');
    sdms{sc} = newsdm;
    newsdm.IncludesConstant = 1;
    newsdm.FirstConfoundPredictor = numel(csn) + 1;
    newsdm.NrOfPredictors = numel(csn) + 1;
    newsdm.NrOfDataPoints = nti;
    newsdm.PredictorColors = floor(255.999 .* rand(numel(csn), 3));
    newsdm.PredictorColors(end+1, :) = 255;
    newsdm.PredictorNames = csn(:)';
    newsdm.PredictorNames{end+1} = 'Constant';
    newsdm.SDMMatrix = xX(sti, csi);
    newsdm.SDMMatrix(:, end + 1) = 1;
    newsdm.RTCMatrix = xX(sti, csi);

    % regressors of interest
    if ~isempty(regoi)

        % match to names
        regii = multimatch(regoi, lower(newsdm.PredictorNames));

        % get those first with a > 0 value in their order
        newidx = regii(regii > 0);
        newcare = numel(newidx);

        % add remaining ones in their order
        newidx = [newidx(:)', setdiff(1:(numel(csn)+1), newidx(:)')];

        % update first confound predictor
        newsdm.FirstConfoundPredictor = newcare + 1;

        % re-order
        newsdm.PredictorColors = newsdm.PredictorColors(newidx, :);
        newsdm.PredictorNames = newsdm.PredictorNames(newidx);
        newsdm.SDMMatrix = newsdm.SDMMatrix(:, newidx);
        newsdm.RTCMatrix = newsdm.RTCMatrix(:, newidx(1:end-1));
    end

    % saving
    try
        if sessno > 1
            newsdm.SaveAs(sprintf('%s/%s_run%d.sdm', sdmp{1}, sdmp{2}, sc));
        else
            newsdm.SaveAs(sdmfile(:)');
        end
    catch ne_eo;
        clearxffobjects(sdms);
        error( ...
            'neuroelf:xffError', ...
            'Error saving SDM to file: ''%s''.', ...
            ne_eo.message ...
        );
    end
end

% clean objects
if nargout < 2 && ...
    sessno > 1
    clearxffobjects(sdms(2:end));
end
sdm = sdms{1};
