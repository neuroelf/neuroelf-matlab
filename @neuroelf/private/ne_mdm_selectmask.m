% FUNCTION ne_mdm_selectmask: select mask file
function ne_mdm_selectmask(varargin)

% Version:  v0.9b
% Build:    11050712
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
ch = ne_gcfg.h.MDM.h;

% pre-set with default
of = ch.MaskFile.String{1};
ch.MaskFile.String{1} = 'Click to select...';

% bring up file selector
try
    msk = {};
    if ~isempty(of) && ...
       ~strcmpi(of, 'click to select...')
        msk = {xff([fileparts(of) '/*.msk'])};
    else
        msk = {xff('*.msk')};
    end
    if isempty(msk{1})
        return;
    end
    if ~isxff(msk{1}, 'msk')
        uiwait(warndlg('Invalid MSK file!', 'NeuroElf - warning', 'modal'));
        clearxffobjects(msk);
        return;
    end
catch ne_eo;
    clearxffobjects(msk);
    uiwait(warndlg(['Error occurred:' char([10, 10]) ne_eo.message], ...
        'NeuroElf - user information', 'modal'));
    return;
end

% get layout and filename, then clear mask object
ml = msk{1}.Layout;
mf = msk{1}.FilenameOnDisk;
msk{1}.ClearObject;

% check in comparison to first functional file
try
    v = {};
    v = {xff(ch.FuncFiles.String{1})};
    if isempty(v{1})
        return;
    end
    if ~isxff(v{1}, 'vtc')
        clearxffobjects(v);
        uiwait(warndlg('Invalid first VTC file!', 'NeuroElf - warning', 'modal'));
        return;
    end
catch ne_eo;
    clearxffobjects(v);
    uiwait(warndlg(['Error occurred:' char([10, 10]) ne_eo.message], ...
        'NeuroElf - user information', 'modal'));
    return;
end

% check layouts
vl = v{1}.Layout;
clearxffobjects(v);
if any(ml([1:3, 5:13]) ~= vl([1:3, 5:13]))
    uiwait(warndlg('MSK and VTC files mismatch.', 'NeuroElf - warning', 'modal'));
    return;
end

% set filename
ch.MaskFile.String{1} = mf;
