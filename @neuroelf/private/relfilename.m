function reltarget = relfilename(sourcefile, targetfile)
% relfilename  - build relative filename from two given files
%
% FORMAT:       rfname = relfilename(sourcefile, targetfile)
%
% Input fields:
%
%       sourcefile  relative filename of the originating file
%       targetfile  relative filename of the target file
%
% Output fields:
%
%       rfname      relative filename between the two
%
% relfilename is useful to build relative filenames for HTML
% compliant links between documents in a folder tree

% Version:  v0.9a
% Build:    10051716
% Date:     May-17 2010, 10:48 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

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
if nargin < 2 || ...
   ~ischar(sourcefile) || ...
   ~ischar(targetfile)
    error( ...
        'neuroelf:BadArgument',...
        'Bad or missing argument.' ...
    );
end

% for empty target file use source dir
if isempty(targetfile)
    [reltarget{1:3}] = fileparts(sourcefile);
    reltarget = [reltarget{2:3}];
    return;
end

% accept empty sourcefile
if isempty(sourcefile)
    sourcefile = 'CURRENT.DIR';
else
    sourcefile = strrep(sourcefile(:)', '\', '/');
end
targetfile = strrep(targetfile(:)', '\', '/');

% get filename parts
schain = splittocell(sourcefile, '/', 1);
tchain = splittocell(targetfile, '/', 1);

% test for ms-dos names
if length(schain) > 1 && ...
    length(tchain) > 1 && ...
   (any(schain{1} == ':') || ...
    any(tchain{1} == ':'))
    if any(tchain{1} == ':')
        if strcmp(schain{1}, tchain{1}) || ...
           (ispc && ...
            strcmpi(schain{1}, tchain{1}))
            schain(1) = [];
            tchain(1) = [];
        else
            reltarget = targetfile;
            return;
        end
    elseif targetfile(1) == '/'
        reltarget = [schain{1} targetfile]; return;
    else
        reltarget = [gluetostring({schain{1:(end-1)}}, '/') '/' ...
                     targetfile];
        return;
    end
else
    if targetfile(1) == '/'
        reltarget = targetfile;
        return;
    elseif sourcefile(1) == '/'
        reltarget = [gluetostring({schain{1:(end-1)}}, '/') '/' ...
                     targetfile];
        return;
    end
end

% remove all leading similarities
while length(schain) > 0 && ...
    length(tchain) > 0 && ...
    (strcmp(schain{1}, tchain{1}) || ...
     (ispc && ...
      strcmpi(schain{1}, tchain{1})))
    schain(1) = [];
    tchain(1) = [];
end

% initialize reltarget
reltarget = '';

% for each remaining folder name go one dir up
while length(schain) > 1
    reltarget = ['../' reltarget];
    schain(1) = [];
end

% add remaining folders and file from target
reltarget = [reltarget gluetostring(tchain, '/')];
