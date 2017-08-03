%NEUROSYNTH_DOWNLOAD  Download NeuroSynth database into PLP object
%   NEUROSYNTH_DOWNLOAD downloads the current database file from GitHub
%   and stores it into the _files/neurosynth/rawdata folder as a PLP object.

% Version:  v1.1
% Build:    16053014
% Date:     May-30 2016, 2:27 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2016, Jochen Weber
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

% create temp folder, and change into it
[tpath, tfolder] = fileparts(tempname);
cd(tpath);
mkdir(tfolder);
cd(tfolder);

% download raw data
fprintf('Downloading NeuroSynth database from GitHub...\n');
!curl https://raw.githubusercontent.com/neurosynth/neurosynth-data/master/current_data.tar.gz > neurosynth_current_data.tar.gz

% unpack raw data
fprintf('Unpacking NeuroSynth database...\n');
!tar -xzf neurosynth_current_data.tar.gz

% create PLP (takes a while...)
n = neuroelf;
fprintf('Creating PLP object... (takes a while!)\n');
plp = n.importplp('database.txt', struct('studycol', 'id', 'tal2x', 'n'));

% load features
fprintf('Loading and parsing features...\n');
f = xff('features.txt', 'ntt');

% get features list
flist = f.Header.pmid(:);
firstf = n.findfirst(~cellfun('isempty', regexpi(flist, '^a')));
flist = flist(firstf:end);

% get PMIDs
pmids = f.Data(:, 1);

% and features table
fvals = sparse(f.Data);
fvals(:, 1:firstf) = [];

% add to PLP
plp.RunTimeVars.Features = flist;
plp.RunTimeVars.FeaturesStudy = pmids;
plp.RunTimeVars.FeaturesValues = fvals;

% and finally, add two featurevalue columns for + (and -) contrast weights
plp.ColumnNames{end+1} = 'FValuePos';
plp.ColumnNames{end+1} = 'FValueNeg';
plp.Points(:, end+1:end+2) = 0;
plp.RunTimeVars.ColumnIsText.FValuePos = false;
plp.RunTimeVars.ColumnIsText.FValueNeg = false;
plp.NrOfColumns = size(plp.Points, 2);

% save
tfile = [neuroelf_path('files') filesep 'neurosynth' filesep ...
    'rawdata' filesep 'neurosynth_' datestr(now, 29) '.plp'];
fprintf('Saving %s...\n', tfile);
plp.SaveAs(tfile);

% remove temp path
cd('..');
system(sprintf('rm -rf %s', tfolder));
cd([neuroelf_path('files') filesep 'neurosynth']);
