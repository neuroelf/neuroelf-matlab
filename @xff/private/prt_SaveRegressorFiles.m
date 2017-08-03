function xo = prt_SaveRegressorFiles(xo, opts)
% PRT::SaveRegressorFiles  - saves onsets into flat text files
%
% FORMAT:       prt.SaveRegressorFiles([opts]);
%
% Input fields:
%
%       filename    output filename (use 'struct' if not to be saved)
%       opts        optional settings
%        .ldelim    line delimiter (default: char(10))
%        .namehead  boolean flag, include name of condition (default: false)
%        .numform   number format (default: '%.3f')
%        .numhead   boolean flag, include number of onsets (default: false)
%        .outpat    output filename pattern (default: '%F_%C.txt')
%        .parsep    force-write additional files for parameters (true)
%        .single    single file for all conditions (enables namehead)
%        .tr        TR (only required if tunits needs conversion)
%        .tspec     time-lengths written as either {'duration'} or 'offset'
%        .tunits    time units, either of 'msec', {'sec'}, 'vol'
%        .volorig0  boolean flag, first volume is 0 (default: true)
%        .weights   store weights as third column (default: true)
%
% No output fields.
%
% Using: asciiwrite, gluetostringc.

% Version:  v1.1
% Build:    16062215
% Date:     Jun-22 2016, 3:16 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2012, 2014, 2016, Jochen Weber
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
asciiwrite    = ne_methods.asciiwrite;
gluetostringc = ne_methods.gluetostringc;

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'prt')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'ldelim') || ~ischar(opts.ldelim) || isempty(opts.ldelim)
    opts.ldelim = char(10);
end
if ~isfield(opts, 'namehead') || ~islogical(opts.namehead) || numel(opts.namehead) ~= 1
    opts.namehead = false;
end
if ~isfield(opts, 'numform') || ~ischar(opts.numform) || isempty(opts.numform) || ...
   ~any(opts.numform(:)' == '%')
    opts.numform = '%.3f';
else
    opts.numform = opts.numform(:)';
end
if ~isfield(opts, 'numhead') || ~islogical(opts.numhead) || numel(opts.numhead) ~= 1
    opts.numhead = false;
end
if ~isfield(opts, 'outpat') || ~ischar(opts.outpat) || isempty(opts.outpat)
    opts.outpat = '%f_%c.txt';
end
opts.outpat = strrep(opts.outpat(:)', '%f', regexprep(xo.F, '\.prt$', '', 'ignorecase'));
if ~isfield(opts, 'parsep') || ~islogical(opts.parsep) || numel(opts.parsep) ~= 1
    opts.parsep = true;
end
if ~isfield(opts, 'single') || ~islogical(opts.single) || numel(opts.single) ~= 1
    opts.single = false;
end
if opts.single && numel(bc.Cond) > 1
    opts.namehead = true;
    opts.numhead = true;
end
if ~isfield(opts, 'tr') || ~isa(opts.tr, 'double') || numel(opts.tr) ~= 1 || ...
    isinf(opts.tr) || isnan(opts.tr) || opts.tr <= 0
    opts.tr = 2;
end
if ~isfield(opts, 'tspec') || ~ischar(opts.tspec) || isempty(opts.tspec) || ...
   ~any('do' == lower(opts.tspec(1)))
    opts.tspec = 'd';
else
    opts.tspec = lower(opts.tspec(1));
end
if ~isfield(opts, 'tunits') || ~ischar(opts.tunits) || isempty(opts.tunits) || ...
   ~any('msv' == lower(opts.tunits(1)))
    opts.tunits = 's';
else
    opts.tunits = lower(opts.tunits(1));
end
if ~isfield(opts, 'volorig0') || ~islogical(opts.volorig0) || numel(opts.volorig0) ~= 1
    opts.volorig0 = true;
end
if ~isfield(opts, 'weights') || ~islogical(opts.weights) || numel(opts.weights) ~= 1
    opts.weights = true;
end
if opts.weights
    cform = [opts.numform '\t' opts.numform '\t' opts.numform '%s'];
else
    cform = [opts.numform '\t' opts.numform '%s'];
end

% time conversion required
if lower(bc.ResolutionOfTime(1)) == 'm'
    if opts.tunits == 's'
        tfac = 0.001;
    elseif opts.tunits == 'v'
        tfac = 1 / opts.tr;
    else
        tfac = 1;
    end
else
    tfac = opts.tr;
end

% get data
conds = bc.Cond;
nconds = numel(conds);

% produce overall array?
if opts.single
    scont = cell(nconds, 1);
end

% generate output
for cc = 1:nconds

    % get condition details
    c = conds(cc);

    % get on/offsets (with global offset)
    oo = tfac .* c.OnOffsets;
    noo = size(oo, 1);

    % convert offsets to duration?
    if opts.tspec == 'd'
        oo(:, 2) = oo(:, 2) - oo(:, 1);
    end

    % get weights
    w = c.Weights;
    if isempty(w) || ~opts.weights
        w = ones(noo, 1);
    elseif opts.parsep
        w = cat(2, ones(noo, 1), w);
    end

    % generate contents
    cplus = double(opts.namehead) + double(opts.numhead);
    ccont = cell(cplus + 1, 1);
    if opts.namehead
        ccont{1} = c.ConditionName{1};
        if opts.numhead
            ccont{2} = sprintf('%d', noo);
        end
    elseif opts.numhead
        ccont{1} = sprintf('%d', noo);
    end

    % iterate across parameters
    for wc = 1:size(w, 2)

        % store weights
        if opts.weights
            cvals = cell(noo, 4);
            for oc = 1:noo
                cvals{oc, 1} = oo(oc, 1);
                cvals{oc, 2} = oo(oc, 2);
                cvals{oc, 3} = w(oc, wc);
            end

            % update condition name?
            if size(w, 2) > 1 && opts.namehead
                ccont{1} = sprintf('%s_x_par%d', c.ConditionName{1}, wc);
            end
        else
            cvals = cell(noo, 3);
            for oc = 1:noo
                cvals{oc, 1} = oo(oc, 1);
                cvals{oc, 2} = oo(oc, 2);
            end
        end
        cvals(:, end) = repmat({opts.ldelim}, noo, 1);
        cvals = cvals';
        ccont{end} = sprintf(repmat(cform, 1, noo), cvals{:});

        % save output
        if ~opts.single
            try
                if opts.namehead
                    asciiwrite(strrep(opts.outpat, '%c', ccont{1}), ...
                        gluetostringc(ccont, opts.ldelim, false));
                elseif opts.parsep
                    if size(w, 2) == 1
                        asciiwrite(strrep(opts.outpat, '%c', c.ConditionName{1}), ...
                            gluetostringc(ccont, opts.ldelim, false));
                    elseif wc == 1
                        asciiwrite(strrep(opts.outpat, '%c', ...
                            sprintf('%s_main', c.ConditionName{1})), ...
                            gluetostringc(ccont, opts.ldelim, false));
                    else
                        asciiwrite(strrep(opts.outpat, '%c', ...
                            sprintf('%s_x_par%d', c.ConditionName{1}, wc - 1)), ...
                            gluetostringc(ccont, opts.ldelim, false));
                    end
                elseif size(w, 2) > 1
                    asciiwrite(strrep(opts.outpat, '%c', ...
                        sprintf('%s_x_par%d', c.ConditionName{1}, wc)), ...
                        gluetostringc(ccont, opts.ldelim, false));
                else
                    asciiwrite(strrep(opts.outpat, '%c', c.ConditionName{1}), ...
                        gluetostringc(ccont, opts.ldelim, false));
                end
            catch xfferror
                rethrow(xfferror);
            end
        else
            if wc == 1
                scont{cc} = gluetostringc(ccont, opts.ldelim, true);
            else
                scont{cc} = [scont{cc}, gluetostringc(ccont, opts.ldelim, true)];
            end
        end
    end
end

% single output
if opts.single
    scont = gluetostringc(scont, opts.ldelim, false);
    try
        asciiwrite(strrep(opts.outpat, '%c', ''), scont);
    catch xfferror
        rethrow(xfferror);
    end
end
