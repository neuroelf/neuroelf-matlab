function xo = aft_Browse(xo, varargin)
% AFT::Browse  - add variable to GUI
%
% FORMAT:       obj.Browse;
%
% No input/output fields.
%
% TYPES: AVA, CMP, DMR, DDT, FMR, FSBF, FSMF, GLM, HEAD, HDR, MGH, MSK, MTC, NLF, SMP, SRF, TOM, TVL, VDW, VMP, VMR, VTC
%
% Note: this function requires GUI being available (figure/uicontrol).

% Version:  v1.1
% Build:    21111013
% Date:     Nov-10 2021, 1:15 PM EST
% Author:   Jochen Weber, NeuroElf.net, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2021, Jochen Weber
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
global ne_gcfg;

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true)
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end

% only for certain files
bc = xo.C;
ft = lower(xo.S.Extensions{1});
if any(strcmp(ft, {'ava', 'cmp', 'dmr', 'fmr', 'fsbf', 'fsmf', 'glm', 'hdr', 'head', ...
    'mgh', 'msk', 'mtc', 'nlf', 'smp', 'srf', 'tom', 'v16', 'vdw', 'vmp', 'vmr', 'vtc'}))

    % simply pass control to GUI
    try
        if isempty(ne_gcfg)
            neuroelf_gui;
        end
        neuroelf_gui('openfile', xo);
    catch xfferror
        rethrow(xfferror);
    end

    % for a sub-set
    if any(strcmp(ft, {'cmp', 'glm', 'vmp'})) && nargin > 1 && ...
       (~strcmp(ft, 'glm') || xo.C.ProjectType < 2) && ...
        isa(varargin{1}, 'double') && ~isempty(varargin{1}) && ...
       ~any(isinf(varargin{1}(:)) | isnan(varargin{1}(:))) && ...
        all(varargin{1}(:) > 0 & varargin{1}(:) == fix(varargin{1}(:)))
        try
            neuroelf_gui('setcstatmap', varargin{1}(:));
        catch xfferror
            fprintf('Error setting current stat map: %s.\n', xfferror.message);
        end
    elseif any(strcmp(ft, {'fsmf', 'glm', 'smp'})) && nargin > 1 && ...
       (~strcmp(ft, 'glm') || xo.C.ProjectType == 2) && ...
        isa(varargin{1}, 'double') && ~isempty(varargin{1}) && ...
       ~any(isinf(varargin{1}(:)) | isnan(varargin{1}(:))) && ...
        all(varargin{1}(:) > 0 & varargin{1}(:) == fix(varargin{1}(:)))
        try
            neuroelf_gui('setcsrfstatmap', varargin{1}(:));
        catch xfferror
            fprintf('Error setting current SRF stat map: %s.\n', xfferror.message);
        end
    elseif strcmp(ft, 'mtc') && nargin > 1 && ...
        isa(varargin{1}, 'double') && ~isempty(varargin{1}) && ...
       ~any(isinf(varargin{1}(:)) | isnan(varargin{1}(:))) && ...
        all(varargin{1}(:) > 0 & varargin{1}(:) == fix(varargin{1}(:)))
        if nargin > 2 && isa(varargin{2}, 'double') && numel(varargin{2}) == 1 && ...
           ~isinf(varargin{2}) && ~isnan(varargin{2}) && varargin{2} >= 1
            bc.RunTimeVars.SubMapVol = varargin{2};
            xo.C = bc;
        end
        try
            neuroelf_gui('setcsrfstatmap', varargin{1}(:));
        catch xfferror
            fprintf('Error setting current MTC map: %s.\n', xfferror.message);
        end
    elseif strcmp(ft, 'vtc') && nargin > 1 && ...
        isa(varargin{1}, 'double') && ~isempty(varargin{1}) && ...
       ~any(isinf(varargin{1}(:)) | isnan(varargin{1}(:))) && ...
        all(varargin{1}(:) > 0 & varargin{1}(:) == fix(varargin{1}(:)))
        try
            if isfield(bc.RunTimeVars, 'AvgVTC') && islogical(bc.RunTimeVars.AvgVTC) && ...
                numel(bc.RunTimeVars.AvgVTC) == 1 && bc.RunTimeVars.AvgVTC
                neuroelf_gui('setcstatmap', varargin{1}(:));
            end
        catch xfferror
            fprintf('Error setting current VTC map: %s.\n', xfferror.message);
        end
    end
end
