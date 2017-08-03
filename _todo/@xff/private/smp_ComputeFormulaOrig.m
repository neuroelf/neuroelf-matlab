function xo = smp_ComputeFormulaOrig(xo, formula, opts)
% SMP::ComputeFormula  - compute a new map from existing SMP maps
%
% FORMAT:       [smp = ] smp.ComputeFormula(formula [, opts])
%
% Input fields:
%
%       formula     string giving a formula, supporting the following
%                   #i -> .Map(i).SMPData
%                   $i -> .Map(opts.mapsel(i)).SMPData
%                   whereas i can be a single number, or a valid range
%                   using the i1:i2 or i1:s:i2 format
%       opts        optional settings
%       .mapsel     sub-selection of maps (for enhanced indexing)
%       .name       set target map name to name
%       .pvalues    flag, if true, convert maps to pvalues
%       .source     map used as template, default first map encountered
%       .target     specify a target map index (otherwise added at end)
%                   - additionally all other Map. subfields are accepted
%
% Output fields:
%
%       smp         SMP with added/replaced map
%
% Using: findfirst, mappvalue, parseformula.

% Version:  v1.1
% Build:    16020212
% Date:     Feb-02 2016, 12:27 PM EST
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

% import neuroelf library
using(neuroelf, 'all');

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'smp') || ...
   ~ischar(formula) || isempty(formula)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
formula = formula(:)';

% get content
bc = xo.C;
maps = bc.Map(:);
nummaps = numel(maps);

% check options
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'mapsel') || ((~isa(opts.mapsel, 'double') || ...
     any(isinf(opts.mapsel(:)) | isnan(opts.mapsel(:)) | opts.mapsel(:) < 1 | ...
         opts.mapsel(:) > nummaps | opts.mapsel(:) ~= fix(opts.mapsel(:)))) && ...
    (~ischar(opts.mapsel) || isempty(opts.mapsel)))
    opts.mapsel = [];
elseif isa(opts.mapsel, 'double')
    opts.mapsel = opts.mapsel(:);
else
    mnames = {maps(:).Name};
    opts.mapsel = ~cellfun('isempty', regexpi(mnames, opts.mapsel(:)'));
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
    isinf(opts.target) || isnan(opts.target) || opts.target < 1 || opts.target > (nummaps + 1)
    opts.target = nummaps + 1;
else
    opts.target = fix(opts.target);
end

% parse formula
fname = 'SMPData';
try
    if any(formula == '$') && ~isempty(opts.mapsel)
        formula = strrep(formula, ':$', sprintf(':%d', numel(opts.mapsel)));
        if ~opts.pvalues
            formula = parseformula(formula, '$', ['bc.Map($).' fname], 4, opts.mapsel);
        else
            formula = parseformula(formula, ...
                '$', ['mappvalue(bc.Map($).' fname ', bc.Map($))'], 4, opts.mapsel);
        end
    end
    if any(formula == '#')
        if ~opts.pvalues
            formula = parseformula(formula, '#', ['bc.Map(#).' fname], 4, 1:nummaps);
        else
            formula = parseformula(formula, ...
                '#', ['mappvalue(bc.Map(#).' fname ', bc.Map(#))'], 4, 1:nummaps);
        end
    end
catch xfferror
    rethrow(xfferror);
end

% valid formula?
fm = strfind(formula, 'bc.Map(');
if isempty(fm)
    error('neuroelf:xff:badArgument', 'Formula didn''t produce any indexing operation.');
end

% get first map if necessary
if isempty(opts.source)
    [fns, fne] = regexp(formula(fm(1):end), '\d+');
    opts.source = str2double(formula((fm(1)-1) + (fns(1):fne(1))));
end

% try to perform computation
try
    newmap = bc.Map(opts.source);
    newmap.([fname 'CT']) = [];
    newmap.RunTimeVars.Formula = oformula;
    newmap.RunTimeVars.FormulaSelection = {opts.mapsel, selmnames, mnames};
    if ~opts.pvalues
        newmap.(fname) = single(eval(formula));
    else
        newmap.(fname) = single(mappvalue(eval(formula), newmap, true));
    end
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
            bc.MapParameter(fc).Values(end+1:opts.target) = ...
                bc.MapParameter(fc).Values(opts.source);
        end
    end
end

% main SMP doesn't have RunTimeVars in Map?
if ~isfield(bc.Map, 'RunTimeVars')
    for fc = 1:numel(bc.Map)
        bc.Map(fc).RunTimeVars = struct;
    end
end

% set into target
bc.Map(opts.target) = newmap;
bc.NrOfMaps = numel(bc.Map);
bc.RunTimeVars.AutoSave = true;
xo.C = bc;
