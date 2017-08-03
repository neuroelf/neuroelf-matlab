function tvmp = onsetshifttest(vtc, prt, cons, range, opts)
% onsetshifttest  - test models with different shifts on VTC
%
% FORMAT:       vmp = onsetshifttest(vtc, prt, cons [, range, [opts]])
%
% Input fields:
%
%       vtc         VTC filename or object
%       prt         PRT filename or object
%       cons        contrast definition(s) (see GLM::FFX_tMap)
%       range       onset shift range in ms (default: -2TR : 1/2TR : 2TR)
%       opts        options used in VTC::CreateGLM
%
% Output fields:
%
%       vmp         VMP object with shifted maps

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:29 AM EST
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

% argument check
if nargin < 3 || ...
   (~isxff(vtc, 'vtc') && ...
    ~ischar(vtc)) || ...
   (~isxff(prt, 'prt') && ...
    ~ischar(prt)) || ...
   ~isa(cons, 'double') || ...
    isempty(cons) || ...
    any(isinf(cons(:)) | isnan(cons(:)))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin < 4 || ...
   ~isa(range, 'double') || ...
    any(isinf(range(:)) | isnan(range(:)))
    range = [];
else
    range = sort(range(:))';
end
if nargin < 5 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end

% load objects as required
if ischar(vtc)
    vtcfile = true;
    try
        vtc = xff(vtc(:)');
        if ~isxff(vtc, 'vtc')
            error( ...
                'neuroelf:InvalidArgument', ...
                'Invalid VTC filename.' ...
            );
        end
    catch ne_eo;
        rethrow(ne_eo);
    end
else
    vtcfile = false;
end
if ischar(prt)
    prtfile = true;
    try
        prt = xff(prt(:)');
        if ~isxff(prt, 'prt')
            error( ...
                'neuroelf:InvalidArgument', ...
                'Invalid PRT filename.' ...
            );
        end
    catch ne_eo;
        rethrow(ne_eo);
    end
else
    prtfile = false;
end

% check range
if isempty(range)
    range = (-2*vtc.TR):(0.5*vtc.TR):(2*vtc.TR);
else
    range(abs(range) > (0.5*vtc.TR*size(vtc.VTCData, 1))) = [];
end

% initiate progress bar
nr = numel(range);
try
    pbar = xprogress;
    xprogress(pbar, 'setposition', [80, 200, 480, 36]);
    xprogress(pbar, 'settitle', 'Shifted model computation...');
    xprogress(pbar, 0, sprintf('%d models...', nr), 'visible', 0, 1);
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    pbar = [];
end

% for each onset in range
for rc = 1:nr

    % make setting
    opts.tshift = range(rc);

    % compute GLM
    glm = vtc.CreateGLM(prt, opts);

    % compute contrast maps
    vmp = glm.FFX_tMap(cons);

    % add shift to contrast names
    for mc = 1:numel(vmp.Map)
        vmp.Map(mc).Name = sprintf('%s (shift: %dms)', vmp.Map(mc).Name, range(rc));
    end

    % keep or update maps
    if rc == 1
        tvmp = vmp;
    else
        tvmp.Map(end+1:end+mc) = vmp.Map;
        vmp.ClearObject;
    end

    % clear GLM object
    glm.ClearObject;

    % progress bar
    if ~isempty(pbar)
        pbar.Progress(rc / nr);
    end
end

% re-order maps
tvmp.Map = tvmp.Map(lsqueeze(reshape(1:(mc*nr), mc, nr)'));

% unload objects
if vtcfile
    vtc.ClearObject;
end
if prtfile
    prt.ClearObject;
end

% clear progress bar
if ~isempty(pbar)
    closebar(pbar);
end
