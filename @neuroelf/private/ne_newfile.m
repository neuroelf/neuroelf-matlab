% FUNCTION ne_newfile: create new object
function varargout = ne_newfile(varargin)
% ne_newfile  - add a new file to workspace and, possibly, controls
%
% FORMAT:       [obj = ] ne_newfile(src, evt, type)
%
% Input fields:
%
%       src         handle of UI object issuing the call (0 for GUI)
%       evt         event data (currently unused)
%       type        either of 'vmp', 'vmr'
%
% Output fields:
%
%       obj         object reference (e.g. for filenames)

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:39 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, 2014, 2016, Jochen Weber
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

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% only with valid input
if nargin < 3 || ...
   ~ischar(varargin{3}) || ...
    numel(varargin{3}) ~= size(varargin{3}, 2) || ...
   ~any(strcmpi(varargin{3}, {'vmp', 'vmr'}))
    return;
end

% create object
try
    newobj = xff(['new:' lower(varargin{3})]);
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    uiwait(warndlg(['Error creating new ' upper(varargin{3}) ': ' ...
        ne_eo.message], 'NeuroElf - error message', 'modal'));
    return;
end
if nargout > 0
    varargout{1} = newobj;
end

% VMR resolution
if lower(varargin{3}(end)) == 'r'

    % get setting
    vres = ne_gcfg.c.ini.MainFig.NewVMRRes;
    if vres < 1
        newobj.VoxResX = vres;
        newobj.VoxResY = vres;
        newobj.VoxResZ = vres;
        fc = 256 * 2 ^ round(-log2(vres));
        newobj.FramingCube = fc;
        newobj.VMRData(fc, fc, fc) = 0;
    end
end

% open file
ne_openfile(0, 0, newobj, true);

% get corresponding variable name
objname = cv_varlist(0, 0, newobj);

% set as tool-loaded
if ~isempty(objname)
    ne_gcfg.wc.(objname) = true;
end
