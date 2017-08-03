function [prt, prts, TR] = spmmat2prt(spmmat, prtfile)
% spmmat2prt  - convert a SPM.mat file into BV's PRT file(s)
%
% FORMAT:       [prt, prts, TR] = spmmat2prt(spmmat, prtfile)
%
% Input fields:
%
%       spmmat      either SPM.mat filename or struct with fields
%        .SPM       1x1 struct containing loaded struct
%             - or -
%        .xBF.UNITS units for timing information
%        .xBF.T     number of time bins (also to calc TR)
%        .xBF.dt    delta time for time bins (event-related designs)
%        .Sess      sessions
%          .U       predictors
%            .name  predictor name
%            .ons   onsets
%            .dur   durations
%            .P     parameters (not yet supported by BV)
%
% Output fields:
%
%       prt         object of first session
%       prts        cell array with other session objects
%       TR          TR (useful for further scripting)

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
   ~ischar(prtfile) || ...
    isempty(prtfile)
    error( ...
        'neuroelf:BadArgument', ...
        'Missing or bad argument. Try ''help %s''.', ...
        mfilename ...
    );
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
   ~isfield(spmmat, 'Sess') || ...
   ~isstruct(spmmat.Sess) || ...
    isempty(spmmat.Sess) || ...
   ~isfield(spmmat.Sess, 'U') || ...
   ~isfield(spmmat, 'xBF') || ...
   ~isstruct(spmmat.xBF) || ...
    isempty(spmmat.xBF) || ...
   ~isfield(spmmat.xBF, 'UNITS') || ...
   ~isfield(spmmat.xBF, 'T') || ...
   ~isfield(spmmat.xBF, 'dt')
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid SPM.mat structure given.' ...
    );
end

% get some settings
bf_t  = spmmat.xBF.T;
bf_dt = spmmat.xBF.dt;
bf_u  = lower(spmmat.xBF.UNITS);
if ~any(strcmp(bf_u(:)', {'scan', 'scans', 'sec', 'secs'}))
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid units found in SPM.mat: %s.', ...
        bf_u(:)' ...
    );
end
TR = bf_t * bf_dt;
sess = spmmat.Sess;
switch (bf_u(:)')
    case {'scan', 'scans'}
        ResolutionOfTime = 'Volumes';
        fac = 1;
    case {'sec', 'secs'}
        ResolutionOfTime = 'msec';
        fac = 1000;
end
sessno = numel(sess);
[prtabs{1:2}] = isabsolute(prtfile(:)');
[prtp{1:3}] = fileparts(prtabs{2});
if isempty(prtp{1})
    prtp{1} = '.';
end

% generate protocols
prts = cell(1, sessno);
for sc = 1:sessno
    newprt = xff('new:prt');
    newprt.FileVersion = 3;
    newprt.ResolutionOfTime = ResolutionOfTime;
    U = sess(sc).U;
    for uc = 1:numel(U)
        uname = U(uc).name;
        if iscell(uname)
            uname = uname{1};
        end
        uons = U(uc).ons * fac;
        udur = U(uc).dur * fac;
        udur(udur == 0) = 1;
        onoff = zeros(numel(uons), 2);

        % volume based
        if fac == 1

            % BrainVoyager is one-based on volumes
            onoff(:, 1) = uons(:) + 1;

        % msec based
        else
            onoff(:, 1) = uons(:);
        end
        onoff(:, 2) = onoff(:, 1) + udur(:) - 1;

        % parametric regression
        if isfield(U(uc), 'P') && ...
            isstruct(U(uc).P) && ...
           ~isempty(U(uc).P) && ...
            isfield(U(uc).P, 'P') && ...
            numel(U(uc).P(1).P) == size(onoff, 1)

            % add to onoff
            upt = 3;
            for upc = 1:numel(U(uc).P)
                try
                    onoff(:, upt) = U(uc).P(upc).P(:);
                    onoff(:, upt) = onoff(:, upt) - mean(onoff(:, upt));
                    upt = upt + 1;
                catch ne_eo;
                    neuroelf_lasterr(ne_eo);
                end
            end
        end

        % add to protocol
        newprt.AddCond(uname, onoff);
    end

    % save PRT
    try
        if sessno > 1
            newprt.Experiment = sprintf('%s - Run %d', prtp{2}, sc);
            newprt.SaveAs(sprintf('%s/%s_run%d.prt', prtp{1}, prtp{2}, sc));
        else
            newprt.Experiment = prtp{2};
            newprt.SaveAs(prtfile(:)');
        end
    catch ne_eo;
        newprt.ClearObject;
        clearxffobjects(prts);
        error( ...
            'neuroelf:xffError', ...
            'Error saving PRT to file: ''%s''.', ...
            ne_eo.message ...
        );
    end

    % only if requested
    if nargout > 1 || ...
        sc == 1
        prts{sc} = newprt;
    else
        newprt.ClearObject;
    end
end

% return first session as prt
prt = prts{1};

% for no output
if nargout < 1
    prt.ClearObject;
end
