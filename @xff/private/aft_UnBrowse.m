function xo = aft_UnBrowse(xo, varargin)
% AFT::UnBrowse  - remove variable from GUI
%
% FORMAT:       obj.UnBrowse;
%
% No input/output fields.
%
% TYPES: DMR, FMR, GLM, HEAD, HDR, MGH, MSK, NLF, SMP, SRF, VDW, VMP, VMR, VTC
%
% Note: this function requires GUI being available (figure/uicontrol).

% Version:  v1.1
% Build:    16031615
% Date:     Mar-16 2016, 3:47 PM EST
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

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true)
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end

% for certain files, do nothing
ft = xo.S.Extensions{1};
if ~any(strcmpi(ft, {'cmp', 'dmr', 'fmr', 'glm', 'head', 'hdr', 'mgh', ...
    'msk', 'nlf', 'smp', 'srf', 'v16', 'vdw', 'vmp', 'vmr', 'vtc'}))
    return;
end

% simply pass control to GUI
try
    neuroelf_gui('closefile', xo);
catch xfferror
    rethrow(xfferror);
end
