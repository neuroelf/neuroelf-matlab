function xo = srf_ApplySPMsn(xo, spmsn)
% SRF::ApplySPMsn  - apply SPM-based _sn.mat file on vertex coordinates
%
% FORMAT:       [srf = ] srf.ApplySPMsn(spmsn);
%
% Input fields:
%
%       srf         folded mesh (e.g. RECOSM mesh)
%       spmsn       SPM-based _sn.mat file (filename or content)
%
% Output fields:
%
%       srf         altered SRF object
%
% Using: applyspmsnc.

% Version:  v1.1
% Build:    16031911
% Date:     Mar-19 2016, 11:29 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
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

% neuroelf library
global ne_methods;

% check arguments
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'srf') || ...
   ((~ischar(spmsn) || isempty(spmsn) || exist(spmsn(:)', 'file') ~= 2) && ...
    (~isstruct(spmsn) || numel(spmsn) ~= 1 || ~isfield(spmsn, 'VG') || ...
     ~isfield(spmsn, 'VF') || ~isfield(spmsn, 'Tr') || ~isfield(spmsn, 'Affine') || ...
     ~isstruct(spmsn.VF) || numel(spmsn.VF) ~= 1 || ~isfield(spmsn.VF, 'mat') || ...
     ~isa(spmsn.VF.mat, 'double') || ~isequal(size(spmsn.VF.mat), [4, 4]) || ...
     ~isstruct(spmsn.VG) || isempty(spmsn.VG) || ~isfield(spmsn.VG, 'dim') || ~isfield(spmsn.VG, 'mat') || ...
     ~isa(spmsn.VG(1).dim, 'double') || ~isequal(size(spmsn.VG(1).dim), [1, 3]) || ...
     ~isa(spmsn.VG(1).mat, 'double') || ~isequal(size(spmsn.VG(1).mat), [4, 4]) || ...
     ~isa(spmsn.Tr, 'double') || ndims(spmsn.Tr) ~= 4 || ...
     ~isa(spmsn.Affine, 'double') || ~isequal(size(spmsn.Affine), [4, 4])))
    error('neuroelf:xff:BadArgument', 'Invalid call to %s.', mfilename);
end

% load SPM sn.mat file
if ischar(spmsn)
    try
        spmsn = load(spmsn);

        % and apply from here
        srf_ApplySPMsn(xo, spmsn);
    catch xfferror
        rethrow(xfferror);
    end

    % then return
    return;
end

% get contents
bc = xo.C;
sh = xo.H;
if isfield(sh, 'SUpdate') && ...
   (~isa(sh.SUpdate, 'function_handle') || numel(sh.SUpdate) ~= 1) && ...
   (~iscell(sh.SUpdate) || isempty(sh.SUpdate) || ~isa(sh.SUpdate{1}, 'function_handle'))
    xo.H = rmfield(xo.H, 'SUpdate');
    sh = xo.H;
end

% apply transform
newcoord = ne_methods.applyspmsnc(128 - bc.VertexCoordinate(:, [3, 1, 2]), spmsn.Tr, ...
    spmsn.VG(1).dim, inv(spmsn.VG(1).mat), spmsn.VF.mat * spmsn.Affine);
bc.VertexCoordinate = 128 - newcoord(:, [2, 3, 1]);

% set back in global memory array
xo.C = bc;

% recalc normals
srf_RecalcNormals(xo);

% update in GUI
if isfield(sh, 'SUpdate')
    if iscell(sh.SUpdate)
        feval(sh.SUpdate{1}, 0, 0, sh.SUpdate{2:end});
    else
        feval(sh.SUpdate);
    end
    drawnow;
end
