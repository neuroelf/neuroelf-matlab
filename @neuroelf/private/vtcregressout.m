function vtcregressout(infile, rfiles, outfile, opts)
% vtcregressout  - regress variance out of a VTC
%
% FORMAT:       vtcregressout([infile [, rfiles [, outfile [, opts]]]])
%
% Input fields:
%
%       infile      input VTC filename
%       rfiles      char or cell array with filename(s)
%       outfile     output filename
%       opts        optional settings
%        .tfiltfrq  number of temporal filtering frequencies (default: 0)
%        .tfilttyp  filtering type, either of 'DCT', {'Fourier'}
%        .trans     perform either 'psc' or 'z' transform (default: 'none')
%
% No output fields.
%
% Note: all non-given arguments will be requested via the GUI!

% Version:  v0.9b
% Build:    11050712
% Date:     Jul-08 2010, 12:52 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, Jochen Weber
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
   ~ischar(infile) || ...
    numel(infile) < 4 || ...
   ~strcmpi(infile(end-3:end), '.vtc')
    infile = '';
end
if nargin < 2 || ...
   (~ischar(rfiles) && ...
    ~iscell(rfiles)) || ...
    isempty(rfiles)
    rfiles = {};
elseif ischar(rfiles)
    rfiles = {rfiles(:)'};
end
for rc = numel(rfiles):-1:1
    if ~ischar(rfiles{rc}) || ...
        numel(rfiles{rc}) < 4
        rfiles(rc) = [];
    end
end
if nargin < 3 || ...
   ~ischar(outfile) || ...
    numel(outfile) < 4 || ...
   ~strcmpi(outfile(end-3:end), '.vtc')
    outfile = '';
end
if nargin < 4 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end

% try to load VTC
vtc = {};
try
    if ~isempty(infile)
        vtc = {xff(infile)};
    else
        vtc = {xff('*.vtc', 'Please select VTC to regress variance out of...')};
    end
    if ~isxff(vtc{1}, 'vtc')
        error( ...
            'neuroelf:BadInput', ...
            'No VTC file selected.' ...
        );
    end
catch ne_eo;
    clearxffobjects(vtc);
    rethrow(ne_eo);
end

% select input files if required
if isempty(rfiles)
    sfiles = true;
    while sfiles
        [sffile, sfpath] = uigetfile( ...
            {'*.sdm', 'BrainVoyager QX single-study design matrix files (*.sdm)'; ...
             '*.rtc', 'BrainVoyager QX reference time course files (*.rtc)'; ...
             '*.txt', 'SPM realignment parameter files (rp*.txt)'; ...
             '*.mat', 'MAT files with one numeric variable (*.mat)'}, ...
            'Please select another file with regressors', 'MultiSelect', 'off');
        if isequal(sffile, 0) || ...
            isequal(sfpath, 0) || ...
           ~ischar(sffile) || ...
            isempty(sffile)
            sfiles = false;
        else
            rfiles{end+1} = [sfpath '/' sffile];
        end
    end
end
if isempty(rfiles)
    warning( ...
        'neuroelf:BadInput', ...
        'No regressor files specified.' ...
    );
end

% regression
vtc{1}.RegressOut(rfiles{:}, opts);

% save as
try
    if isempty(outfile)
        vtc{1}.SaveAs;
    else
        vtc{1}.SaveAs(outfile(:)');
    end
catch ne_eo;
    clearxffobjects(vtc);
    error( ...
        'neuroelf:SavingFailed', ...
        'Error saving VTC file: %s.', ...
         ne_eo.message ...
    );
end
clearxffobjects(vtc);
