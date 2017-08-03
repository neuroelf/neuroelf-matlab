function xo = aft_ComputeFormula(xo, formula, opts)
% AFT::ComputeFormula  - compute a new map from existing data
%
% FORMAT:       [obj = ] obj.ComputeFormula(formula [, opts])
%
% Input fields:
%
%       formula     string giving a formula, supporting the following
%                   @#i.Field -> {.Map(i).Field}
%                   @$i.Field -> {.Map(opts.mapsel(i)).Field}
%                   #i -> e.g. .Map(i).VMPData for VMPs
%                   $i -> e.g. .Map(opts.mapsel(i)).VMPData for VMPs
%                   whereas i can be a single number, or a valid range
%                   using the i1:i2 or i1:s:i2 format
%       opts        optional settings
%       .mapsel     sub-selection of maps (for enhanced indexing using $i)
%       .name       set target map name to name
%       .pvalues    flag, if true, convert maps to pvalues
%       .source     map used as template, default first map encountered
%       .target     specify a target map index (otherwise added at end)
%                   - additionally all other Map. subfields are accepted
%
% Output fields:
%
%       vmp         VMP with added/replaced map
%
% TYPES: HDR, SMP, VMP, VTC.
%
% Using: findfirst, lsqueeze, mappvalue, multimatch, parseformula.

% Version:  v1.1
% Build:    16021321
% Date:     Feb-13 2016, 9:29 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010 - 2015, 2016, Jochen Weber
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

% import functions from neuroelf library (allow ALL functions to be used!)
using(neuroelf, 'all');

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, {'hdr', 'smp', 'vmp', 'vtc'}) || ...
   ~ischar(formula) || isempty(formula)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
formula = formula(:)';

% get content
ftype = lower(xo.S.Extensions{1});
bc = xo.C;
if ftype(end) == 'p'
    maps = bc.Map(:);
    fname = [upper(ftype) 'Data'];
elseif ftype(end) == 'r'
    maps = bc.RunTimeVars.Map(:);
    fname = 'VoxelData';
else
    maps = repmat(struct('Type', 15, 'Name', 'TimePoint', 'DF1', 1, ...
        'LowerThreshold', 0, 'UpperThreshold', double(max(bc.VTCData(:)))), ...
        size(bc.VTCData, 1), 1);
    fname = 'VTCData';
end
nummaps = numel(maps);

% check options
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
mnames = lsqueeze({maps(:).Name});
if ~isfield(opts, 'mapsel') || ((~isa(opts.mapsel, 'double') || ...
    any(isinf(opts.mapsel(:)) | isnan(opts.mapsel(:)) | opts.mapsel(:) < 1 | ...
         opts.mapsel(:) > nummaps | opts.mapsel(:) ~= fix(opts.mapsel(:)))) && ...
    (~ischar(opts.mapsel) || isempty(opts.mapsel)))
    opts.mapsel = [];
elseif isa(opts.mapsel, 'double')
    opts.mapsel = opts.mapsel(:);
else
    opts.mapsel = ~cellfun('isempty', regexpi(mnames, opts.mapsel(:)'));
end
if ~isempty(opts.mapsel)
    selmnames = mnames(opts.mapsel);
else
    selmnames = {};
end
if ~isfield(opts, 'name') || ~ischar(opts.name) || numel(opts.name) ~= size(opts.name, 2)
    opts.name = formula;
end
if ~isfield(opts, 'pvalues') || ~islogical(opts.pvalues) || numel(opts.pvalues) ~= 1
    opts.pvalues = false;
end
if ~isfield(opts, 'source') || ~isa(opts.source, 'double') || numel(opts.source) ~= 1 || ...
    isinf(opts.source) || isnan(opts.source) || opts.source < 1 || opts.source > nummaps || ...
    opts.source ~= fix(opts.source)
    opts.source = [];
end
if ~isfield(opts, 'target') || ~isa(opts.target, 'double') || numel(opts.target) ~= 1 || ...
    isinf(opts.target) || isnan(opts.target) || opts.target < 0 || opts.target > (nummaps + 1)
    opts.target = nummaps + 1;
else
    opts.target = fix(opts.target);
end

% parse formula
try
    oformula = formula;
    formula = strrep(formula, ':$', sprintf(':%d', numel(opts.mapsel)));
    formula = strrep(formula, ':#', sprintf(':%d', numel(opts.mapsel)));
    formula = regexprep(formula, '@\$(\d+\:\d+\:\d+)\.([A-Za-z][_A-Za-z0-9]*)', ...
        '{maps(opts.mapsel($1)).$2}');
    formula = regexprep(formula, '@\$(\d+\:\d+)\.([A-Za-z][_A-Za-z0-9]*)', ...
        '{maps(opts.mapsel($1)).$2}');
    formula = regexprep(formula, '@\$(\d+)\.([A-Za-z][_A-Za-z0-9]*)', ...
        'maps(opts.mapsel($1)).$2');
    formula = regexprep(formula, '#\$(\d+\:\d+\:\d+)\.([A-Za-z][_A-Za-z0-9]*)', ...
        '{maps($1).$2}');
    formula = regexprep(formula, '#\$(\d+\:\d+)\.([A-Za-z][_A-Za-z0-9]*)', ...
        '{maps($1).$2}');
    formula = regexprep(formula, '#\$(\d+)\.([A-Za-z][_A-Za-z0-9]*)', ...
        'maps($1).$2');
    if any(formula == '$') && ~isempty(opts.mapsel)
        if ~opts.pvalues
            if ftype(end) == 'p'
                formula = parseformula(formula, '$', ['bc.Map($).' fname], 4, opts.mapsel);
            elseif ftype(end) == 'r'
                formula = parseformula(formula, ...
                    '$', ['bc.' fname '(:, :, :, $)'], 4, opts.mapsel);
            else
                formula = parseformula(formula, ...
                    '$', ['permute(bc.' fname '($, :, :, :), [2, 3, 4, 1])'], ...
                    4, opts.mapsel);
            end
        else
            if ftype(end) == 'p'
                formula = parseformula(formula, ...
                    '$', ['mappvalue(bc.Map($).' fname ', maps($))'], 4, opts.mapsel);
            elseif ftype(end) == 'r'
                formula = parseformula(formula, ...
                    '$', ['mappvalue(bc.' fname '(:, :, :, $), maps($))'], 4, opts.mapsel);
            else
                formula = parseformula(formula, ...
                    '$', ['mappvalue(permute(bc.' fname '($, :, :, :), [2, 3, 4, 1]), maps($))'], ...
                    4, opts.mapsel);
            end
        end
    end
    if any(formula == '#')
        if ~opts.pvalues
            if ftype(end) == 'p'
                formula = parseformula(formula, '#', ['bc.Map(#).' fname], 4, 1:nummaps);
            elseif ftype(end) == 'r'
                formula = parseformula(formula, ...
                    '#', ['bc.' fname '(:, :, :, #)'], 4, 1:nummaps);
            else
                formula = parseformula(formula, ...
                    '#', ['permute(bc.' fname '(#, :, :, :), [2, 3, 4, 1])'], 4, 1:nummaps);
            end
        else
            if ftype(end) == 'p'
                formula = parseformula(formula, ...
                    '#', ['mappvalue(bc.Map(#).' fname ', maps(#))'], 4, 1:nummaps);
            elseif ftype(end) == 'r'
                formula = parseformula(formula, ...
                    '#', ['mappvalue(bc.' fname '(:, :, :, #), maps(#))'], 4, 1:nummaps);
            else
                formula = parseformula(formula, ...
                    '#', ['mappvalue(permute(bc.' fname '(#, :, :, :), [2, 3, 4, 1]), maps(#))'], 4, 1:nummaps);
            end
        end
    end
catch xfferror
    rethrow(xfferror);
end

% valid formula?
fm = strfind(formula, 'bc.');
if isempty(fm)
    error('neuroelf:xff:badArgument', 'Formula didn''t produce any indexing operation.');
end

% get first map if necessary
if isempty(opts.source)
    [fns, fne] = regexp(formula(fm(1):end), '\d+');
    opts.source = str2double(formula((fm(1)-1) + (fns(1):fne(1))));
end

% new object
if opts.target < 1
    xo = aft_CopyObject(xo);
    opts.target = 1;
    if ftype(end) == 'p'
        bc.Map = bc.Map(1);
    elseif ftype(end) == 'r'
        bc.VoxelData = bc.VoxelData(:, :, :, 1);
        bc.VoxelDataCT = {[]};
    else
        bc.VTCData = bc.VTCData(1, :, :, :);
    end
elseif ftype(end) == 'r' && istransio(bc.VoxelData)
    bc.VoxelData = resolve(bc.VoxelData);
elseif ftype(end) == 'c' && istransio(bc.VTCData)
    bc.VTCData = resolve(bc.VTCData);
end

% try to perform computation
try
    if ftype(end) == 'p'
        newmap = bc.Map(opts.source);
        newmap.([fname 'CT']) = [];
    elseif ftype(end) == 'r'
        newmap = bc.RunTimeVars.Map(opts.source);
    else
        newmap = maps(opts.source);
    end
    if ~opts.pvalues
        newmapdata = single(eval(formula));
    else
        newmapdata = single(mappvalue(eval(formula), newmap, true));
    end
    if ftype(end) == 'p'
        newmap.(fname) = newmapdata;
    end
    newmap.Name = opts.name;
    newmap.RunTimeVars.Formula = oformula;
    newmap.RunTimeVars.FormulaSelection = {opts.mapsel, selmnames, mnames};
catch xfferror
    error('neuroelf:xff:badArgument', 'Error computing formula: ''%s''.', xfferror.message);
end

% set fields
of = fieldnames(opts);
mf = fieldnames(newmap);
for fc = 1:numel(mf)
    fm = findfirst(strcmpi(of, mf{fc}));
    if ~isempty(fm)
        newmap.(mf{fc}) = opts.(of{fm});
    end
end

% also extend MapParameter if needed
if isfield(bc, 'MapParameter') && ~isempty(bc.MapParameter)
    for fc = 1:numel(bc.MapParameter)
        if numel(bc.MapParameter(fc).Values) < opts.target
            bc.MapParameter(fc).Values(end+1:opts.target) = bc.MapParameter(fc).Values(opts.source);
        end
    end
end

% main ?MP doesn't have RunTimeVars in Map?
if ftype(end) == 'p' && ~isfield(bc.Map, 'RunTimeVars')
    for fc = 1:numel(bc.Map)
        bc.Map(fc).RunTimeVars = struct;
    end
end

% extend bc.RunTimeVars.Map fieldnames otherwise (as necessary)
if ftype(end) ~= 'p'
    ofnames = fieldnames(bc.RunTimeVars.Map);
    nfnames = fieldnames(newmap);
    if numel(ofnames) ~= numel(nfnames) || ~all(strcmp(ofnames, nfnames))
        nfnames = nfnames(multimatch(nfnames, ofnames) == 0);
        for fc = 1:numel(nfnames)
            bc.RunTimeVars.Map(1).(nfnames{fc}) = [];
        end
    end
end

% set into target
if ftype(end) == 'p'
    bc.Map(opts.target) = newmap;
    bc.NrOfMaps = numel(bc.Map);
elseif ftype(end) == 'r'
    bc.VoxelData(:, :, :, opts.target) = newmapdata;
    bc.ImgDim.Dim(5) = size(bc.VoxelData, 4);
    bc.RunTimeVars.Map(opts.target) = newmap;
else
    bc.VTCData(opts.target, :, :, :) = reshape(newmapdata, [1, size(newmapdata)]);
    bc.NrOfVolumes = size(bc.VTCData, 1);
    bc.RunTimeVars.Map(opts.target) = newmap;
end
bc.RunTimeVars.AutoSave = true;
xo.C = bc;
