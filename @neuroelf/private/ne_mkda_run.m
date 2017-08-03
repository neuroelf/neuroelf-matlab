% PUBLIC FUNCTION ne_mkda_run: run MKDA analysis
function varargout = ne_mkda_run(varargin)

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:38 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, 2016, Jochen Weber
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
hFig = ne_gcfg.h.MKDA.MKDAFig;
ch = ne_gcfg.h.MKDA.h;

% only allow one instance
if any(strcmp(ne_gcfg.c.blockcb, 'mkda_run'))
    return;
end
ne_gcfg.c.blockcb{end+1} = 'mkda_run';

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% get content of PLP dropdown
plps = ch.PLPs;
plpud = plps.UserData;
plpid = plps.Value;
try
    plp = plpud{plpid, 3};
    if numel(plp) ~= 1 || ...
       ~isxff(plp, 'plp')
        error( ...
            'neuroelf:GUI:BadPLP', ...
            'Bad PLP object.' ...
        );
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'mkda_run')) = [];
    return;
end
rtv = plp.RunTimeVars;
anaidx = ch.Analyses.Value;
anas = plp.MKDAAnalyses;
mcfg = rtv.Config;
x = xff;

% hide dialog
hFig.Visible = 'off';

% run all
if nargin > 2 && ...
    ischar(varargin{3}) && ...
    strcmpi(varargin{3}(:)', 'runall')

    % iterate over analyses
    sumvmp = [];
    for ac = 1:size(anas, 1)

        % set value
        ch.Analyses.Value = ac;

        % set analysis
        ne_mkda_setana;

        % try to run
        ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'mkda_run')) = [];
        try
            [mvmp, svmp] = ne_mkda_run(0, 0, 'nobrowse');

        % error occurred
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;

            % show error
            uiwait(errordlg(ne_eo.message, 'NeuroElf - error', 'modal'));

            % show figure
            hFig.Visible = 'on';
            figure(hFig.MLHandle);

            % browse partial map
            if numel(sumvmp) == 1 && ...
                isxff(sumvmp, 'vmp')
                ne_openfile(0, 0, sumvmp, true);
            end

            % return
            return;
        end

        % browse individual maps
        ne_openfile(0, 0, mvmp, true);

        % copy maps over
        if isempty(sumvmp)
            if ~isempty(svmp)
                sumvmp = svmp;
            end
        else
            sumvmp.Map = catstruct(sumvmp.Map(:)', svmp.Map(:)');
            svmp.ClearObject;
        end
    end

    % browse map
    if ~isempty(sumvmp)
        ne_openfile(0, 0, sumvmp, true);
    end

% run the selected analysis
else

    % get analysis configuration
    ana = anas{anaidx, 2};

    % create configuration
    mkdacfg = struct( ...
        'applymask', mcfg.ApplyMask, ...
        'asimiter',  ana.Iterations, ...
        'asimkeep',  true, ...
        'asimkthr',  true, ...
        'asimmask',  [], ...
        'asimsmpl',  mcfg.SpatialNull, ...
        'asimrthr',  ne_gcfg.c.ini.Tools.alphasim.Thresholds(:)', ...
        'bbox',      [], ...
        'cond',      ch.Points.UserData{1}, ...
        'contcomp',  mcfg.ContrastComp, ...
        'contexclw', mcfg.ContrastCompExclWeight, ...
        'contnames', {{ana.Contrast}}, ...
        'contnull',  ana.NullDist, ...
        'contrasts', {{ana.Contrast}}, ...
        'grpmeth',   mcfg.GroupMapComp, ...
        'indivmaps', mcfg.KeepIndivMaps, ...
        'jbmeth',    mcfg.JoinBlobComp, ...
        'pbar',      [], ...
        'res',       ana.Resolution, ...
        'scale',     ana.Scaling, ...
        'smkern',    ana.SphereSize, ...
        'smkinterp', 'linear', ...
        'smkmdist',  ana.SphereTaper, ...
        'smkres',    1, ...
        'studycol',  ana.StudyColumn, ...
        'stwf',      ana.Weights, ...
        'stwp',      mcfg.PPSWeighting, ...
        'stwsel',    true, ...
        'unique',    mcfg.UniqueUnitPoints, ...
        'usecons',   findfirst(strcmpi(ana.ContColumn, plp.ColumnNames(:))), ...
        'usesize',   [], ...
        'usevalue',  []);

    % masking
    msko = [];
    mskl = false;
    try
        msko = x.Document(ana.Mask);
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
        if any(ana.Mask == '/')
            try
                msko = xff(ana.Mask);
                mskl = true;
            catch ne_eo;
                ne_gcfg.c.lasterr = ne_eo;
            end
        else
            try
                msko = neuroelf_file('c', ana.Mask);
                mskl = true;
            catch ne_eo;
                ne_gcfg.c.lasterr = ne_eo;
            end
        end
    end
    mkdacfg.asimmask = msko;

    % run MKDA
    try
        mvmp = plp.MKDA(mkdacfg);
        if mskl
            msko.ClearObject;
        end
    catch ne_eo;
        if mskl
            clearxffobjects({msko});
        end
        ne_gcfg.c.lasterr = ne_eo;
        if nargin < 3 || ...
           ~ischar(varargin{3}) || ...
           ~strcmpi(varargin{3}(:)', 'nobrowse')
            uiwait(errordlg(ne_eo.message, 'NeuroElf - error', 'modal'));
            hFig.Visible = 'on';
            figure(hFig.MLHandle);
            ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'mkda_run')) = [];
            return;
        else
            ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'mkda_run')) = [];
            rethrow(ne_eo);
        end
    end

    % set as output
    if nargout > 0
        varargout{1} = mvmp;
    end

    % copy summary VMPs
    svmp = [];
    if mcfg.SummaryVMP && ...
        mcfg.KeepIndivMaps
        svmp = mvmp.CopyObject;
        smaps = svmp.Map;
        smaptype = cat(1, smaps.Type);
        svmp.Map(smaptype == 16) = [];
        if nargout > 1
            varargout{2} = svmp;
        end
    elseif ~mcfg.KeepIndivMaps
        smaps = mvmp.Map;
        smaptype = cat(1, smaps.Type);
        mvmp.Map(smaptype == 16) = [];
    end

    % single run of multiple runs?
    if nargin > 2 && ...
        ischar(varargin{3}) && ...
        strcmpi(varargin{3}(:)', 'nobrowse')

        % simply return output after unblocking next run
        ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'mkda_run')) = [];
        return;
    end

    % otherwise, open files
    ne_openfile(0, 0, mvmp, true);
    if mcfg.SummaryVMP && ...
        isxff(svmp, 'vmp')
        ne_openfile(0, 0, svmp, true);
    end
end

% set the dialog visible again
hFig.Visible = 'on';
figure(hFig.MLHandle);

% re-enable next analysis run
ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'mkda_run')) = [];
