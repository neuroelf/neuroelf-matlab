% FUNCTION ne_vmp_applysvc: apply SVC to a VMP map
function varargout = ne_vmp_applysvc(varargin)

% Version:  v0.9d
% Build:    14071116
% Date:     Jul-11 2014, 4:39 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2014, Jochen Weber
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

% only valid if single VMP and map
stvar = ne_gcfg.fcfg.StatsVar;
stvix = ne_gcfg.fcfg.StatsVarIdx;
if numel(stvar) ~= 1 || ...
   ~isxff(stvar, 'vmp') || ...
    numel(stvix) ~= 1
    return;
end

% what SVC source
if nargin > 2 && ...
    ischar(varargin{3}) && ...
    strcmpi(varargin{3}(:)', 'voi')
    if isempty(ne_gcfg.h.Clusters.Value)
        return;
    end
    svc = ne_gcfg.voi;
    if numel(svc) ~= 1 || ...
        ~isxff(svc, 'voi') || ...
        isempty(svc.VOI)
        return;
    end
    svc = svc.CopyObject;
    svc.VOI = svc.VOI(ne_gcfg.h.Clusters.Value);
    svc.NrOfVOIs = numel(svc.VOI);
elseif nargin > 2 && ...
    ischar(varargin{3}) && ...
    strcmpi(varargin{3}(:)', 'vmr')
    if ~isxff(ne_gcfg.fcfg.SliceVar, {'hdr', 'head', 'msk', 'vmr', 'vtc'})
        return;
    end
    svc = ne_gcfg.fcfg.SliceVar.CopyObject;
elseif nargin > 2 && ...
    ischar(varargin{3}) && ...
   ~isempty(varargin{3}) && ...
    exist(varargin{3}(:)', 'file') == 2
    try
        svc = [];
        svc = xff(varargin{3}(:)');
        if ~isxff(svc, {'hdr', 'head', 'msk', 'voi', 'vmr', 'vtc'})
            error( ...
                'neuroelf:BadObject', ...
                'Invalid object for SVC.' ...
            );
        end
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
        clearxffobjects({svc});
        return;
    end
else
    [svcfile, svcpath] = uigetfile({ ...
        '*.msk', 'BrainVoyager Mask files (*.msk)'; ...
        '*.hdr;*.nii;*.nii.gz', 'Analyze/NIftI files (*.hdr, *.nii, *.nii.gz)'; ...
        '*.head', 'AFNI HEAD (+BRIK) files (*.head)'; ...
        '*.voi', 'BrainVoyager VOI files (*.voi)'}, 'Please select SVC mask...');
    if isequal(svcfile, 0) || ...
        isequal(svcpath, 0) || ...
        isempty(svcfile)
        return;
    end
    if isempty(svcpath)
        svcpath = pwd;
    end
    svcfile = [svcpath '/' svcfile];
    if exist(svcfile, 'file') ~= 2
        return;
    end
    try
        svc = [];
        svc = xff(svcfile);
        if ~isxff(svc, {'hdr', 'head', 'msk', 'voi'});
            error( ...
                'neuroelf:BadObject', ...
                'Invalid object for SVC.' ...
            );
        end
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
        clearxffobjects({svc});
    end
end

% copy map
stvar.Map(end+1) = stvar.Map(stvix);
newmap = numel(stvar.Map);

% apply SVC
mfp = ne_gcfg.h.MainFig.Pointer;
try
    ne_gcfg.h.MainFig.Pointer = 'watch';
    drawnow;
    stvar.ClusterTable(newmap, [], struct('svc', svc, 'minsize', 1));
catch ne_eo;
    uiwait(warndlg(ne_eo.message, 'NeuroElf - error', 'modal'));
end
ne_gcfg.h.MainFig.Pointer = mfp;

% delete SVC object
clearxffobjects({svc});

% update
stvar.Browse(newmap);
drawnow;
