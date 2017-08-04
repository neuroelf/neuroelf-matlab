function [b, bnl, v] = neuroelf_build
% neuroelf_build  - return the build number of NeuroElf
%
% FORMAT:       [build, bnl, v] = neuroelf_build
%
% No input fields.
%
% Output fields:
%
%       b           build number (year/month/day/hour)
%       bnl         Fx3 cell array file filenames, build, and version
%       v           version of NeuroElf (from neuroelf_build)
%
% Using: findfiles, grep.

% Version:  v1.1
% Build:    17080412
% Date:     Aug-04 2017, 12:48 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2016, 2017, Jochen Weber
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

% no argument check necessary

% neuroelf object and function handles
nelf = neuroelf;
findfiles = nelf.findfiles;
grep = nelf.grep;

% get main build
[b, v] = getbuildversion(which('neuroelf_build'), grep);

% if more outputs required
if nargout > 1
    p = neuroelf_path;
    toolboxfiles = [ ...
        findfiles(p, '*.m',   'relative='); ...
        findfiles(p, '*.h',   'relative='); ...
        findfiles(p, '*.c',   'relative='); ...
        findfiles(p, '*.cpp', 'relative='); ...
        findfiles(p, '*.*ff', 'relative='); ...
        findfiles(p, '*.tfg', 'relative=')];
    bnl = cell(numel(toolboxfiles), 3);
    for fc = 1:size(bnl, 1)
        bnl{fc, 1} = toolboxfiles{fc};
        [bnl{fc, 2:3}] = getbuildversion([p filesep toolboxfiles{fc}], grep);
    end
end

% sub function to extract build information:
function [b, v] = getbuildversion(filename, grep)
try
    b = str2double(regexprep(grep(filename, ['Bui' 'ld']), ...
        ['^.*Bui' 'ld\:\s*(\d+)\s*.*$'], '$1'));
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    b = NaN;
    v = '?.?';
    return;
end
if nargout > 1
    v = regexprep(grep(filename, 'Version'), ...
        '^.*Version\:\s*v(\d+\.\d+[a-z]*)\s*.*$', '$1');
    if iscell(v)
        if ~isempty(v)
            v = v{1};
        else
            v = '';
        end
    end
end
