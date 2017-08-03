function [t, tfile] = analyzerkfile
% analyzerkfile  - analyze a recorded-keys file
%
% FORMAT:       [t, tfile] = analyzerkfile
%
% No input fields.
%
% Output fields:
%
%       t           time spent on one of three states
%       tfile       selected filename

% Version:  v0.9d
% Build:    14062015
% Date:     Jun-20 2014, 3:51 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2012 - 2014, Jochen Weber
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

% output defaults to 0
t = zeros(2, 3);
tfile = '';

% get filename
[file, folder] = uigetfile({'*.mat','Record-Keys file (*.mat)'}, ...
    'Please select a MAT-file with recordkeys output...');
if isequal(file, 0) || ...
    isequal(folder, 0)
    return;
end

% try to load file
if isempty(folder)
    folder = pwd;
end
tfile = [folder '/' file];
try
    file = load(tfile);
catch ne_eo;
    rethrow(ne_eo);
end

% correct field in file?
if ~isfield(file, 'id') || ...
   ~iscell(file.id) || ...
    isempty(file.id) || ...
    ndims(file.id) ~= 2 || ...
    size(file.id, 2) ~= 2 || ...
   ~all(cellfun(@iscell, file.id(:, 1)))
    return;
end
for cc = 1:size(file.id, 1)
    if numel(file.id{cc, 1}) ~= 1 || ...
       ~ischar(file.id{cc, 1}{1})
        error('analyserkfile:error', 'Bad id cell %d.', cc);
    end
    file.id{cc, 1} = file.id{cc, 1}{1};
end

% find begin and end markers
np0 = find(strcmp(file.id(:, 1), 'numpad0'));
if numel(np0) < 2
    error('No start and end marker found');
end

% get all numpads inbetween (assuming first item is a 3!)
nps = file.id(np0(1):np0(2)-1, 1);
nps{1} = 'numpad3';
nps = strrep(nps, 'numpad', '');

% and corresponding times (difference only)
npt = diff(cat(1, file.id{np0(1):np0(2), 2}));

% which task
validtask = (strcmp(nps, '1') | strcmp(nps, '2') | strcmp(nps, '3'));
if ~all(validtask)
    warning('analyserkfile:warning', 'Unknown key pressed');
    nps(~validtask) = [];
    npt(~validtask) = [];
end

% convert
nps(strcmp(nps, '1')) = {1};
nps(strcmp(nps, '2')) = {2};
nps(strcmp(nps, '3')) = {3};
nps = cat(1, nps{:});

% total time
tt = sum(npt);

% time spent in each task
t1 = sum(npt(nps == 1));
t2 = sum(npt(nps == 2));
t3 = sum(npt(nps == 3));

% output
t = [t1, t2, t3; ([t1, t2, t3] ./ tt)];
