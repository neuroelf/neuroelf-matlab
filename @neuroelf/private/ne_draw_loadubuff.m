% FUNCTION ne_draw_loadubuff: (re-)load UndoBuffer of SliceVar
function ne_draw_loadubuff(varargin)

% Version:  v0.9d
% Build:    14060816
% Date:     Jun-08 2014, 4:00 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2014, Jochen Weber
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
cc = ne_gcfg.fcfg;

% check SliceVar
svar = cc.SliceVar;
if numel(svar) ~= 1 || ...
   ~isxff(svar, {'hdr', 'vmr'})
    return;
end
svartyp = lower(svar.Filetype);

% load VMR for UndoBuffer
udbv = xff('*.(vmr|hdr)', 'Please select VMR/Nifti file containing UndoBuffer for current object...');
if numel(udbv) ~= 1 || ...
   ~isxff(udbv, true)
    return;
end
if ~isxff(udbv, svartyp)
    udbv.ClearObject;
    return;
end

% for VMRs
if strcmp(svartyp, 'vmr')

    % check datatype
    if ~strcmp(class(udbv.VMRData(1)), class(svar.VMRData(1)))
        udbv.ClearObject;
        uiwait(warndlg('VMRs inconsistent in datatype.', 'NeuroElf - info', 'modal'));
        return;
    end

    % check resolution
    if ~all([svar.VoxResX, svar.VoxResY, svar.VoxResZ] == ...
            [udbv.VoxResX, udbv.VoxResY, udbv.VoxResZ])
        udbv.ClearObject;
        uiwait(warndlg('VMRs inconsistent in resolution.', 'NeuroElf - info', 'modal'));
        return;
    end

    % get onsets and size
    tsz = size(svar.VMRData);
    [o1, o2, sz] = bboverlay(udbv, svar);

    % no overlap
    if any(sz < 1)
        udbv.ClearObject;
        uiwait(warndlg('VMRs do not overlap.', 'NeuroElf - info', 'modal'));
        return;
    end

    % ensure UndoBuffer is valid
    if ~isfield(svar.RunTimeVars, 'UndoBuffer') || ...
       ~isequal(size(svar.RunTimeVars.UndoBuffer), tsz)
        svar.RunTimeVars.UndoBuffer = svar.VMRData(1);
        svar.RunTimeVars.UndoBuffer(1) = [];
        svar.RunTimeVars.UndoBuffer(tsz(1), tsz(2), tsz(3)) = 0;
    end

    % copy part we need
    svar.RunTimeVars.UndoBuffer( ...
        o2(1)+1:o2(1)+sz(1), o2(2)+1:o2(2)+sz(2), o2(3)+1:o2(3)+sz(3)) = ...
        udbv.VMRData(o1(1)+1:o1(1)+sz(1), o1(2)+1:o1(2)+sz(2), o1(3)+1:o1(3)+sz(3));

% for HDR
else
    uiwait(warndlg('Not yet implemented.', 'NeuroElf - info', 'modal'));
end

% clear vmr
udbv.ClearObject;
