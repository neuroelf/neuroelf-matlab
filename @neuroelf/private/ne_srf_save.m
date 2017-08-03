function varargout = ne_srf_save(varargin)
% ne_srf_save  - save currently selected SurfVar
%
% FORMAT:       [srf = ] ne_srf_save(SRC, EVT [, filename])
%
% Input fields:
%
%       SRC, EVT    Matlab handle callback inputs (discarded)
%       filename    target filename if to be saved under different name
%
% Output fields:
%
%       srf         object handle of saved SRF
%
% Notes:
%
%     instead of using a filename, the special token 'saveas' can be
%     passed into the function to call AFT::SaveAs for the current
%     surface object.
%
% Examples:
%
%     ne_srf_save(0, 0, 'subject123_LH_Inflated.srf');

% Version:  v1.1
% Build:    16031521
% Date:     Mar-15 2016, 9:53 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2014, 2016, Jochen Weber
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
ch = ne_gcfg.h;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% get/check current SRF
srf = cc.SurfVar;
if numel(srf) ~= 1 || ...
   ~isxff(srf, {'fsbf', 'srf'})
    return;
end
if nargout > 0
    varargout{1} = srf;
end

% no further arguments
e = [];
if nargin < 3 || ...
   ~ischar(varargin{3}) || ...
    numel(varargin{3}) < 4

    % echo
    if ne_gcfg.c.echo
        ne_echo('srf', 'Save');
    end

    % try saving to disk
    try
        if ~isempty(srf.FilenameOnDisk)
            srf.Save;
        else
            srf.SaveAs;
        end
    catch ne_eo;
        e = ne_eo;
    end

% filename given
elseif any(strcmpi(varargin{3}(end-3:end), {'fsbf', '.srf'}))

    % echo
    if ne_gcfg.c.echo
        ne_echo('srf', sprintf('SaveAs(%s)', varargin{3}(:)'))
    end

    % save
    try
        srf.SaveAs(varargin{3}(:)');
    catch ne_eo;
        e = ne_eo;
    end

% general "saveas" call
elseif strcmpi(varargin{3}(:)', 'saveas')

    % echo
    if ne_gcfg.c.echo
        ne_echo('srf', 'SaveAs');
    end

    % try saving as
    try
        srf.SaveAs;
    catch ne_eo;
        e = ne_eo;
    end
end

% error handing
if ~isempty(e)
    uiwait(warndlg(e.message, 'NeuroElf - warning', 'modal'));

% otherwise
else

    % update filename
    sidx = ch.SurfVar.Value;
    svars = ch.SurfVar.UserData;
    [fp, fn, fe] = fileparts(srf.FilenameOnDisk);
    fn = [fn, fe];
    if isempty(fe)
        return;
    end
    srfh = handles(srf);
    if sidx <= size(svars, 1) && ...
        isxff(svars{sidx, 4}, {'fsbf', 'srf'}) && ...
        srf == svars{sidx, 4}
        svarnames = ch.SurfVar.String;
        if ischar(svarnames)
            svarnames = cellstr(svarnames);
        end
        svarnames{sidx} = fn;
        ch.SurfVar.String = svarnames;
        ch.SurfVar.UserData{sidx, 1} = fn;
    end
    svars = ch.Scenery.UserData;
    for sidx = 1:size(svars, 1)
        if isxff(svars{sidx, 4}, {'fsbf', 'srf'}) && ...
            srf == svars{sidx, 4}
            svarnames = ch.Scenery.String;
            if ischar(svarnames)
                svarnames = cellstr(svarnames);
            end
            if isempty(srfh.VertexMorphMeshes)
                svarnames{sidx} = sprintf('%s (NrOfVertices: %d)', fn, srf.NrOfVertices);
            else
                svarnames{sidx} = sprintf( ...
                    '%s (NrOfVertices: %d, NrOfMorphStates: %d)', fn, ...
                    srf.NrOfVertices, size(srfh.VertexMorphMeshes, 1) + 1);
            end
            svars{sidx, 1} = svarnames{sidx};
            ch.Scenery.String = svarnames;
            ch.Scenery.UserData = svars;
            break;
        end
    end
end
