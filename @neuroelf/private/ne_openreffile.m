% FUNCTION ne_openreffile: open referenced file
function ne_openreffile(varargin)

% Version:  v0.9b
% Build:    11112514
% Date:     Apr-10 2011, 4:52 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, Jochen Weber
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

% global variable
global ne_gcfg;
ch = ne_gcfg.h;
try
    cu = ch.(varargin{3});
catch ne_eo;
    ne_gcfg.c.laster = ne_eo;
    return;
end
ud = cu.UserData;

% close any currently open file
if ~isempty(ud{1})
    try
        ne_closefile(0, 0, ud{1});
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
    ud{1} = [];
    cu.UserData = ud;
end

% get currently selected file
csel = cu.Value;

% open the file
if csel > 1
    try
        ud{1} = ne_openfile(0, 0, ud{csel});
        if ~isxff(ud{1}, true)
            ud{1} = [];
        elseif isxff(ud{end}, 'glm')
            glm = ud{end};
            st = glm.Study(csel - 1);
            if isfield(st, 'RunTimeVars') && ...
                isstruct(st.RunTimeVars)
                ud{1}.SetHandle('StudyData', st.RunTimeVars);
                ne_setslicepos;
            end
        end

        % and put reference into UserData for later
        cu.UserData = ud;
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
end
