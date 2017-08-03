% PUBLIC FUNCTION ne_savecluster: save current cluster information
function varargout = ne_savecluster(varargin)

% Version:  v1.1
% Build:    16040916
% Date:     Apr-09 2016, 4:23 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2016, Jochen Weber
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

% save as MSK?
if nargin > 2 && ischar(varargin{3}) && strcmpi(varargin{3}(:)', 'asmask')
    voi = ne_gcfg.voi;
    if ~isxff(voi, 'voi') || isempty(voi.VOI)
        return;
    end
    vvoi = voi.VOI;
    vidx = ne_gcfg.h.Clusters.Value;
    if isempty(vidx)
        vidx = 1:numel(vvoi);
    end
    vcrd = cat(1, vvoi.Voxels);
    if isempty(vcrd)
        return;
    end
    try
        msk = [];
        msk = xff('new:msk');
        stvar = ne_gcfg.fcfg.StatsVar;
        if isxff(stvar, {'ava', 'cmp', 'glm', 'vmp'})
            msk.Resolution = stvar.Resolution;
            msk.XStart = stvar.XStart;
            msk.YStart = stvar.YStart;
            msk.ZStart = stvar.ZStart;
            msk.XEnd = stvar.XEnd;
            msk.YEnd = stvar.YEnd;
            msk.ZEnd = stvar.ZEnd;
            msk.Mask = uint8(zeros(stvar.GetVolumeSize));
        elseif isxff(ne_gcfg.fcfg.SliceVar, {'msk', 'vtc'})
            slvar = ne_gcfg.fcfg.SliceVar;
            msk.Resolution = slvar.Resolution;
            msk.XStart = slvar.XStart;
            msk.YStart = slvar.YStart;
            msk.ZStart = slvar.ZStart;
            msk.XEnd = slvar.XEnd;
            msk.YEnd = slvar.YEnd;
            msk.ZEnd = slvar.ZEnd;
            msk.Mask = uint8(zeros(slvar.GetVolumeSize));
        else
            if size(vcrd, 1) > numel(vidx)
                resdet = abs(diff(vcrd));
                resdet = resdet(:);
                resdet(resdet == 0) = [];
                resdet = min(resdet);
                while resdet > 3
                    resdet = round(0.5 * resdet);
                end
            else
                resdet = 1;
            end
            msk.Resolution = resdet;
            msk.XStart = 44;
            msk.YStart = 38;
            msk.ZStart = 44;
            msk.XEnd = 242;
            msk.YEnd = 194;
            msk.ZEnd = 212;
            msk.Mask = uint8(zeros(round([198, 156, 168] ./ resdet)));
        end
        msk.Mask(bvcoordconv(vcrd, 'tal2bvx', msk.BoundingBox)) = 1;
        msk.SaveAs;
        ne_openfile(0, 0, msk, true);
        return;
    catch ne_eo;
        if isxff(msk, true)
            msk.ClearObject;
        end
        uiwait(warndlg(ne_eo.message, 'NeuroElf - error', 'modal'));
        return;
    end
end

% let xff do the job
try
    ne_gcfg.voi.SaveAs;
    if ne_gcfg.c.echo
        ne_echo('voi.SaveAs;');
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
end
