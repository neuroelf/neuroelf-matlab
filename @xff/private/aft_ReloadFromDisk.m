function xo = aft_ReloadFromDisk(xo)
% AFT::ReloadFromDisk  - try reloading the object from disk
%
% FORMAT:       obj.ReloadFromDisk;
%
% No input / output fields.
%
% TYPES: ALL
%
% Note: if the object requires large amount of memory, this method
%       is bound to fail!

% Version:  v1.1
% Build:    16020214
% Date:     Feb-02 2016, 2:29 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2014, 2016, Jochen Weber
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

% global config
global xffsngl;

% only valid for single file
if numel(xo) ~= 1 || ~xffisobject(xo, true) || xo.L(1) == 'X'
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end

% must have a filename
if isempty(xo.F) || ~isabsolute(xo.F) || exist(xo.F, 'file') ~= 2
    error('neuroelf:xff:badFileName', 'Cannot reload from disk, invalid file: ''%s''.', xo.F);
end

% try to load file
rs = xffsngl.CONF.reloadsame;
try
    nobj = [];
    xffsngl.CONF.reloadsame = true;
    nobj = xff(xo.F, xo.S.Extensions{1});
    if ~xffisobject(nobj, true, xo.S.Extensions{1})
        error('neuroelf:xff:badContent', 'Invalid file content.');
    end
catch xfferror
    xffsngl.CONF.reloadsame = rs;
    delete(nobj);
    rethrow(xfferror);
end
xffsngl.CONF.reloadsame = rs;

% get object contents then clear re-loaded object
objcont = nobj.C;
delete(nobj);

% keep RunTimeVars as they were!
objcont.RunTimeVars = xo.C.RunTimeVars;

% copy contents from reloaded object to calling object
xo.C = objcont;
