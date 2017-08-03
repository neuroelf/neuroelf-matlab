function map = aft_RSA(xo, opts)
% AFT::RSA  - compute a RSA (representational similarity analysis)
%
% FORMAT:       out = aft.RSA(opts)
%
% Input fields:
%
%       opts        structure with mandatory fields
%        .cells     1xC cell arrays with 1xM map/contrast names for cells
%                   and optional fields
%        .slrad     search-light radius (default: 2 * res)
%        .subsel    subject selection (numeric list)
%        .voi       VOI object with regions to test (instead of searchlight)
%
% Output fields:
%
%       map         VMP object with C maps
%
% Using: calcbetas, conjval, findfirst, fitrobustbisquare_img, glmtstat,
%        limitrangec, lsqueeze, resestsmooth, robustnsamplet_img, robustt,
%        sdist, varc, ztrans.

% Version:  v1.1
% Build:    16020214
% Date:     Feb-02 2016, 2:20 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2014, 2016, Jochen Weber
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
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, {'glm', 'vmp'}) || ...
   ~isstruct(opts) || numel(opts) ~= 1 || ~isfield(opts, 'cells') || ~iscell(opts.cells) || ...
    numel(opts.cells) < 2 || isempty(opts.cells{1}) || isempty(opts.cells{2})
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end

% test cells
ncells = opts.cells;
for cc = 1:ncells
    if isempty(opts.cells{cc}) || ~iscell(opts.cells{cc}) || numel(opts.cells{cc}) < 2 || ...
       ~all(cellfun(@ischar, opts.cellc{cc}(:))) || any(cellfun('isempty', opts.cells(:)))
        error('neuroelf:xff:badArgument', 'Invalid cell specified.');
    end
end

% get file content
bc = xo.C;
glmfile = xo.F;
glmid = xo.L;
if isempty(glmfile)
    glmfile = glmid;
end

% VOI
voigiven = false;
if isfield(opts, 'voi') && numel(opts.voi) == 1 && xffisobject(opts.voi, true, 'voi')
    voic = opts.voi.C;
    if ~isempty(voic.VOI)
        voigiven = true;
    end
end

% for GLMs
if strcmpi(xo.Extensions{1}, 'glm')
    if bc.ProjectTypeRFX == 0
        error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
    end
    isrfx = (bc.ProjectTypeRFX > 0);
    ffxspred = glm_SubjectPredictors(xo);
    if isrfx
        numsubs = bc.NrOfSubjects;
        numspred = bc.NrOfSubjectPredictors;
    else
        ffxpred = bc.Predictor;
        ffxpred = {ffxpred(:).Name2};
        ffxpred = ffxpred(:);
        ffxsubs = glm_Subjects(xo);
        numsubs = numel(ffxsubs);
        numspred = numel(ffxspred) + 1;
    end
    if ~isfield(opts, 'subsel') || ~isa(opts.subsel, 'double') || ...
        isempty(opts.subsel) || any(isinf(opts.subsel(:)) | isnan(opts.subsel(:)))
        ga = ones(numsubs, 1);
    else
        opts.subsel = opts.subsel(:)';
        opts.subsel(opts.subsel < 1 | opts.subsel > numsubs) = [];
        ga = zeros(numsubs, 1);
        ga(unique(round(opts.subsel))) = 1;
    end
    if isrfx
        szmap = size(bc.GLMData.RFXGlobalMap);
    else
        szmap = size(bc.GLMData.MCorrSS);
    end

% for VMPs
else

end

% generate output -> VMP (searchlight)
sdist = ne_methods.sdist;
if ~voigiven
    map = xff('new:vmp');
    mapc = map.C;
    mapc.XStart = bc.XStart;
    mapc.XEnd = bc.XEnd;
    mapc.YStart = bc.YStart;
    mapc.YEnd = bc.YEnd;
    mapc.ZStart = bc.ZStart;
    mapc.ZEnd = bc.ZEnd;
    mapc.Resolution = bc.Resolution;
    mapc.RunTimeVars.TrfPlus = bc.RunTimeVars.TrfPlus;
    mapc.Map.Type = 1;
    mapc.Map.LowerThreshold = -sdist('tinv', 0.005, nval);
    mapc.Map.UpperThreshold = -sdist('tinv', 0.0001, nval);
    mapc.Map.Name = 'Contrast:';
    mapc.Map.NrOfLags = [];
    mapc.Map.MinLag = [];
    mapc.Map.MaxLag = [];
    mapc.Map.CCOverlay = [];
    mapc.Map.ClusterSize = 25;
    mapc.Map.EnableClusterCheck = 0;
    mapc.Map.DF1 = nval;
    mapc.Map.DF2 = 0;
    mapc.Map.BonferroniValue = bc.NrOfVoxelsForBonfCorrection;
    mapc.Map.VMPData = single(zeros(szmap));
    mapc.Map(1:(nummaps*gamx*nmf)) = mapc.Map(1);
end

%         % iterate over contrasts
%         occ = 1;
%         findfirst = ne_methods.findfirst;
%         for cc = 1:size(c, 2)
%
%             % get contrast
%             conc = c(:, cc);
%             coni = find(conc ~= 0);
%
%             % initialize temp map
%             tmpmp = zeros([szmap, numsubs]);
%
%             % allow for subjects with missing data (FFX)
%             keepsubs = true(numsubs, 1);
%             gaxk = gax;
%
%             % fill contrast
%             for pc = coni(:)'
%                 for sc = 1:numsubs
%                     if isrfx
%                         tmpmp(:, :, :, sc) = tmpmp(:, :, :, sc) + conc(pc) .* ...
%                             bc.GLMData.Subject(gax(sc)).BetaMaps(:, :, :, pc);
%                     else
%                         keepsubi = findfirst(~cellfun('isempty', regexpi(ffxpred, ...
%                             sprintf('^subject\\s+%s:\\s*%s', ...
%                             ffxsubs{gax(sc)}, ffxspred{pc}))));
%                         if ~isempty(keepsubi)
%                             tmpmp(:, :, :, sc) = tmpmp(:, :, :, sc) + conc(pc) .* ...
%                                 bc.GLMData.BetaMaps(:, :, :, keepsubi);
%                         else
%                             keepsubs(sc) = false;
%                         end
%                     end
%                 end
%             end
%
%             % remove bad subjects from list
%             keepsubs = keepsubs & ~ne_methods.lsqueeze(all(all(all(isnan(tmpmp), 3), 2), 1));
%             if ~all(keepsubs)
%                 tmpmp(:, :, :, ~keepsubs) = [];
%                 gaxk(~keepsubs) = [];
%             end
%             goodsubs = sum(keepsubs);
%             if ngrp < 1
%                 nval = numel(gaxk) - 1;
%             else
%                 nval = numel(gaxk) - ngrp;
%             end
%             if nval < 1
%                 continue;
%             end
%             mapc.Map(occ).DF1 = nval;
%
%             % set additional data
%             artv = struct( ...
%                 'SourceGLM',   glmfile, ...
%                 'SourceGLMID', glmid, ...
%                 'Contrast',    conc, ...
%                 'Covariates',  opts.covs, ...
%                 'Groups',      {opts.groups}, ...
%                 'MeanRem',     opts.meanr, ...
%                 'RFXGLM',      isrfx, ...
%                 'Robust',      opts.robust, ...
%                 'SubPreds',    {ffxspred}, ...
%                 'SubSel',      gaxk);
%
%             % generate 2nd-level single beta-t-maps
%             if ngrp < 1
%                 if ~opts.robust
%                     if numcovs == 0
%                         mmap = mean(tmpmp, 4);
%                         tmap = sqrt(goodsubs) * (mmap ./ std(tmpmp, 0, 4));
%                         if opts.estfwhm
%                             ptc = repmat(mmap, [1, 1, 1, size(tmpmp, 4)]);
%                         end
%                     else
%                         x = [ones(goodsubs, 1), opts.covs(keepsubs, :)];
%                         [bmaps, ixx, ptc, se] = ne_methods.calcbetas(x, tmpmp, 4);
%                         tmap = ne_methods.glmtstat([1, zeros(1, numcovs)], bmaps, ixx, se);
%                     end
%                     if opts.estfwhm
%                         [artv.FWHMResEst, artv.FWHMResImg] = ...
%                             ne_methods.resestsmooth(tmpmp - ptc, bc.Resolution);
%                     end
%                 else
%                     x = [ones(goodsubs, 1), opts.covs(keepsubs, :)];
%                     [bmaps, wmaps] = ne_methods.fitrobustbisquare_img(x, tmpmp);
%                     tmap = ne_methods.robustt(x, tmpmp, bmaps, wmaps, [1, zeros(1, numcovs)]);
%                     if opts.estfwhm
%                         ptc = zeros(size(tmpmp));
%                         for bmc = 1:size(x, 2)
%                             ptc = ptc + repmat(bmaps(:, :, :, bmc), [1, 1, 1, size(x, 1)]) .* ...
%                                 repmat(reshape(x(:, bmc), [1, 1, 1, size(x, 1)]), szmap);
%                         end
%                         ptc = wmaps .* ptc + (1 - wmaps) .* tmpmp;
%                         [artv.FWHMResEst, artv.FWHMResImg] = ...
%                             ne_methods.resestsmooth(tmpmp - ptc, bc.Resolution);
%                     end
%                 end
%                 tmap(isinf(tmap) | isnan(tmap)) = 0;
%
%                 % set name and map data
%                 mapc.Map(occ).Name = sprintf('%s%s', opts.names{cc}, tmapr);
%                 mapc.Map(occ).VMPData = single(tmap);
%                 mapc.Map(occ).RunTimeVars = artv;
%                 occ = occ + 1;
%
%                 % create additional map without residual
%                 if opts.meanr
%
%                     % compute std==1-scaled data
%                     resmp = tmpmp;
%                     bmaps = 1 ./ sqrt(ne_methods.varc(tmpmp, 4));
%                     bmaps(isinf(bmaps) | isnan(bmaps) | bmaps > 1) = 1;
%                     for sc = 1:goodsubs
%                         resmp(:, :, :, sc) = resmp(:, :, :, sc) .* bmaps;
%                     end
%                     resmp(isinf(resmp) | isnan(resmp)) = 0;
%
%                     % average according to tmap
%                     tmin = -sdist('tinv', 0.25, nval);
%                     tmax = -sdist('tinv', 0.05, nval);
%                     wmp = ne_methods.limitrangec(tmax - (1 / (tmax - tmin)) .* tmap, 0.25, 1, 0);
%                     wmp = wmp .* wmp;
%                     if isequal(size(wmp), size(opts.meanrmsk))
%                         wmp = wmp .* opts.meanrmsk;
%                     end
%                     wmp(tmap == 0) = 0;
%                     for sc = 1:goodsubs
%                         resmp(:, :, :, sc) = resmp(:, :, :, sc) .* wmp;
%                     end
%                     resmp = ne_methods.lsqueeze(sum(sum(sum(resmp, 1), 2), 3)) ./ sum(wmp(:));
%
%                     % non-robust remodeling
%                     x = [ones(goodsubs, 1), ne_methods.ztrans(resmp)];
%                     if ~opts.robust
%
%                         % use calcbetas and glmtstat
%                         [bmaps, ixx, ptc, se] = ne_methods.calcbetas(x, tmpmp, 4);
%                         tmap = ne_methods.glmtstat([1, 0], bmaps, ixx, se);
%                         if opts.estfwhm
%                             [artv.FWHMResEst, artv.FWHMResImg] = ...
%                                 ne_methods.resestsmooth(tmpmp - ptc, bc.Resolution);
%                         end
%
%                     % robust remodeling
%                     else
%                         [bmaps, wmaps] = ne_methods.fitrobustbisquare_img(x, tmpmp);
%                         tmap = ne_methods.robustt(x, tmpmp, bmaps, wmaps, [1, 0]);
%                         if opts.estfwhm
%                             ptc = zeros(size(tmpmp));
%                             for bmc = 1:size(x, 2)
%                                 ptc = ptc + repmat(bmaps(:, :, :, bmc), [1, 1, 1, size(x, 1)]) .* ...
%                                     repmat(reshape(x(:, bmc), [1, 1, 1, size(x, 1)]), szmap);
%                             end
%                             ptc = wmaps .* ptc + (1 - wmaps) .* tmpmp;
%                             [artv.FWHMResEst, artv.FWHMResImg] = ...
%                                 ne_methods.resestsmooth(tmpmp - ptc, bc.Resolution);
%                         end
%                     end
%                     tmap(isinf(tmap) | isnan(tmap)) = 0;
%
%                     % set name and map data
%                     mapc.Map(occ).Name = sprintf('%s%s (without average res.)', ...
%                         opts.names{cc}, tmapr);
%                     mapc.Map(occ).VMPData = single(tmap);
%                     mapc.Map(occ).DF1 = nval - 1;
%                     mapc.Map(occ).LowerThreshold = -sdist('tinv', 0.005, nval - 1);
%                     mapc.Map(occ).UpperThreshold = -sdist('tinv', 0.0001, nval - 1);
%                     mapc.Map(occ).RunTimeVars = artv;
%                     occ = occ + 1;
%
%                     % compute contrast for actual average residual
%                     if ~opts.robust
%                         tmap = ne_methods.glmtstat([0, 1], bmaps, ixx, se);
%                     else
%                         tmap = ne_methods.robustt(x, tmpmp, bmaps, wmaps, [0, 1]);
%                     end
%                     tmap(isinf(tmap) | isnan(tmap)) = 0;
%
%                     % set name and map data
%                     mapc.Map(occ).Name = sprintf('%s%s (corr. average res.)', ...
%                         opts.names{cc}, tmapr);
%                     mapc.Map(occ).VMPData = single(tmap);
%                     mapc.Map(occ).DF1 = nval - 1;
%                     mapc.Map(occ).LowerThreshold = -sdist('tinv', 0.005, nval - 1);
%                     mapc.Map(occ).UpperThreshold = -sdist('tinv', 0.0001, nval - 1);
%                     mapc.Map(occ).RunTimeVars = artv;
%                     occ = occ + 1;
%
%                     % compute contrast vs. correlation with residual map
%                     if ~opts.robust
%                         tmap = ne_methods.conjval(ne_methods.glmtstat([1, 0.5], bmaps, ixx, se), ...
%                             ne_methods.glmtstat([1, -0.5], bmaps, ixx, se));
%                     else
%                         tmap = ne_methods.conjval(ne_methods.robustt(x, tmpmp, bmaps, wmaps, [1, 0.5]), ...
%                             ne_methods.robustt(x, tmpmp, bmaps, wmaps, [1, -0.5]));
%                     end
%                     tmap(isinf(tmap) | isnan(tmap)) = 0;
%
%                     % set name and map data
%                     mapc.Map(occ).Name = sprintf('%s%s (w/o+xmsk average res.)', ...
%                         opts.names{cc}, tmapr);
%                     mapc.Map(occ).VMPData = single(tmap);
%                     mapc.Map(occ).DF1 = nval - 1;
%                     mapc.Map(occ).LowerThreshold = -sdist('tinv', 0.005, nval - 1);
%                     mapc.Map(occ).UpperThreshold = -sdist('tinv', 0.0001, nval - 1);
%                     mapc.Map(occ).RunTimeVars = artv;
%                     occ = occ + 1;
%                 end
%
%             % for multiple groups
%             else
%                 if ~opts.robust
%                     for gc1 = 1:ngrp
%                         ga1 = find(ga(gaxk) == gc1);
%                         m1 = (1 / numel(ga1)) .* sum(tmpmp(:, :, :, ga1), 4);
%                         s1 = (1 / numel(ga1)) .* var(tmpmp(:, :, :, ga1), [], 4);
%                         for gc2 = (gc1 + 1):ngrp
%                             ga2 = find(ga(gaxk) == gc2);
%                             m2 = (1 / numel(ga2)) .* sum(tmpmp(:, :, :, ga2), 4);
%                             s2 = (1 / numel(ga2)) .* var(tmpmp(:, :, :, ga2), [], 4);
%                             t12 = (m1 - m2) ./ sqrt(s1 + s2);
%                             df12 = ((s1 + s2) .^ 2) ./ ...
%                                 ((1 / (numel(ga1) - 1)) .* s1 .* s1 + ...
%                                  (1 / (numel(ga2) - 1)) .* s2 .* s2);
%                             badv = find(isinf(t12) | isnan(t12) | isnan(df12) | df12 < 1);
%                             t12(badv) = 0;
%                             df12(badv) = 1;
%                             rdf12 = sum(ga == gc1 | ga == gc2) - 2;
%                             mapc.Map(occ).Name = sprintf('%s (%s > %s)', ...
%                                 opts.names{cc}, opts.groups{gc1, 1}, ...
%                                 opts.groups{gc2, 1});
%                             mapc.Map(occ).VMPData = single(...
%                                 sdist('tinv', sdist('tcdf', t12, df12), rdf12));
%                             mapc.Map(occ).DF1 = rdf12;
%                             mapc.Map(occ).RunTimeVars = artv;
%                             mapc.Map(occ).RunTimeVars.SubSel(gas ~= gc1 & gas ~= gc2) = [];
%                             occ = occ + 1;
%                         end
%                     end
%                 else
%                     tmap = ne_methods.robustnsamplet_img(tmpmp, ga(gaxk));
%                     for gmc = 1:gamx
%                         [g2, g1] = find(gag == gmc);
%                         mapc.Map(occ).Name = sprintf('%s%s (%s > %s)', ...
%                             opts.names{cc}, tmapr, opts.groups{g1, 1}, ...
%                             opts.groups{g2, 1});
%                         mapc.Map(occ).VMPData = single(tmap(:, :, :, gmc));
%                         mapc.Map(occ).DF1 = sum(ga == g1 | ga == g2) - 2;
%                         mapc.Map(occ).RunTimeVars = artv;
%                         mapc.Map(occ).RunTimeVars.SubSel(gas ~= g1 & gas ~= g2) = [];
%                         occ = occ + 1;
%                     end
%                 end
%             end
%         end
%         if occ <= numel(mapc.Map)
%             mapc.Map(occ:end) = [];
%         end
%
%     % MTCs
%     case {2}
%
%         % SRF required
%         numvert = bc.NrOfVertices;
%         if ipo
%             if ~isfield(opts, 'srf') || ...
%                 numel(opts.srf) ~= 1 || ...
%                ~xffisobject(opts.srf, true, 'srf')
%                 error( ...
%                     'xff:BadArgument', ...
%                     'Missing or bad SRF reference in map options.' ...
%                 );
%             end
%             srfc = opts.srf.C;
%             if size(srfc.Neighbors, 1) ~= numvert
%                 error( ...
%                     'xff:BadArgument', ...
%                     'Number of vertices mismatch.' ...
%                 );
%             end
%             nei = srfc.Neighbors(:, 2);
%             nnei = numel(nei);
%             neil = zeros(nnei, 12);
%             neis = zeros(nnei, 12);
%             for nc = 1:nnei
%                 neis(nc) = numel(nei{nc});
%                 neil(nc, 1:neis(nc)) = nei{nc};
%             end
%             neii = (neil > 0);
%         end
%
%         map = xff('new:smp');
%         mapc = map.C;
%         mapc.FileVersion = 4;
%         mapc.NrOfVertices = numvert;
%         if ipo
%             mapc.NameOfOriginalSRF = srfs.F;
%         end
%         mapc.Map.Type = 1;
%         mapc.Map.NrOfLags = [];
%         mapc.Map.MinLag = [];
%         mapc.Map.MaxLag = [];
%         mapc.Map.CCOverlay = [];
%         mapc.Map.ClusterSize = 25;
%         mapc.Map.EnableClusterCheck = 0;
%         mapc.Map.DF1 = nval;
%         mapc.Map.DF2 = 0;
%         mapc.Map.BonferroniValue = bc.NrOfVoxelsForBonfCorrection;
%         mapc.Map.Name = 'Contrast:';
%         mapc.Map.SMPData = single(zeros(numvert, 1));
%         mapc.Map(1:nummaps) = mapc.Map(1);
%
%         % iterate over contrasts
%         for cc = 1:size(c, 2)
%
%             % get contrast
%             conc = c(:, cc);
%             coni = find(conc ~= 0);
%
%             % initialize temp map
%             tmpmp = zeros(numvert, numsubs);
%             for pc = coni(:)'
%                 for sc = 1:numsubs
%                     tmpmp(:, sc) = tmpmp(:, sc) + conc(pc) * ...
%                         bc.GLMData.Subject(sc).BetaMaps(:, pc);
%                 end
%             end
%
%             % generate 2nd-level single beta-t-maps
%             if ~opts.robust
%                 tmap = sqrt(numsubs) * (mean(tmpmp, 2) ./ std(tmpmp, 0, 2));
%             else
%                 [bmaps, wmaps] = ne_methods.fitrobustbisquare_img(ones(numsubs, 1), tmpmp);
%                 tmap = ne_methods.robustt(ones(numsubs, 1), tmpmp, bmaps, wmaps, 1);
%             end
%             tmap(isinf(tmap) | isnan(tmap)) = 0;
%
%             % interpolate ?
%             if ipo
%                 rmap = zeros(numvert, 1);
%                 for nc = 1:size(neii, 2)
%                     rmap(neii(:, nc)) = rmap(neii(:, nc)) + ...
%                         tmap(neil(neii(:, nc), nc));
%                 end
%                 tmap = rmap ./ sum(neii, 2);
%             end
%
%             % set name and map data
%             mapc.Map(cc).Name = sprintf('%s%s', opts.names{cc}, tmapr);
%             mapc.Map(cc).SMPData = tmap;
%         end
%
%      % unknown types
%     otherwise
%         error( ...
%             'xff:UnknownOption', ...
%             'Unknown GLM ProjectType: %d', ...
%             bc.ProjectType ...
%         );
% end
%
% % put back
% mapc.NrOfMaps = numel(mapc.Map);
map.C = mapc;
