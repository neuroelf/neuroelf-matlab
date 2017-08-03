function display(S)
% neuroelf::display  - display info on class object

% Version:  v0.9d
% Build:    14082013
% Date:     Aug-20 2014, 1:40 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2014, Jochen Weber
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

% display methods
neos = struct(S);
disp('neuroelf - methods object; supported methods are');
disp(' ');
meths = fieldnames(neos.meth);
for mc = 1:numel(meths)
    mname = meths{mc};
    mins = neos.meth.(mname){5};
    if ~isempty(mins)
        mins = sprintf('%s, ', mins{:});
        mins = ['(' mins(1:end-2) ')'];
    else
        mins = '';
    end
    mouts = neos.meth.(mname){3};
    if ~isempty(mouts)
        if numel(mouts) > 1 || ...
            strcmp(mouts{end}, 'varargout')
            mouts = sprintf('%s, ', mouts{:});
            mouts = ['[' mouts(1:end-2) '] = '];
        else
            mouts = [mouts{1} ' = '];
        end
    else
        mouts = '';
    end
    fprintf('%32s%s%s\n', mouts, mname, mins);
end
