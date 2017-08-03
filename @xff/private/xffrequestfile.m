function fnames = xffrequestfile(nrf, pattern, ext, formats, saveas)
% xff::_requestfile  - request a file interactively
%
% THIS FUNCTION IS AN INTERNAL FUNCTION OF THE CLASS
%
% @xff
%
% AND IT SHOULD NEVER BE CALLED FROM OUTSIDE THE CLASS CODE

% Version:  v1.1
% Build:    16021216
% Date:     Feb-12 2016, 4:31 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
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

% persistent array with extensions to present for '*' extension
persistent xrf_ext;
if isempty(xrf_ext)
    xrf_ext = struct;
    ftl = { ...
        'vmr', 'vmp', ...
        'glm', 'vtc', 'msk', ...
        'fmr', 'amr', 'map', 'dmr', ...
        'hdr', 'head', 'nii', ...
        'srf', 'mtc', 'smp'};
    xrf_ext.l = [true, false(1, numel(ftl) + 1)];
    xrf_ext.n = neuroelf;
    xrf_ext.x = cell(numel(ftl) + 2, 2);
    ftll = xrf_ext.n.gluetostring(ftl, ';*.');
    xrf_ext.x(1, :) = ...
        {['*.' ftll], ['All NeuroElf files (*.' strrep(ftll, ';', ', ') ')']};
    for tc = 1:numel(ftl)
        fmt = ext.(ftl{tc});
        spec = formats.(lower(fmt{1}(end-2:end)))(fmt{2});
        exf = find(strcmpi(spec.Extensions, ftl{tc}));
        desc = spec.Description{exf(1)};
        if isempty(strfind(desc, '*.'))
            xrf_ext.x(tc + 1, :) = {['*.' ftl{tc}], [desc ' (*.' ftl{tc} ')']};
        else
            [exm{1:3}] = regexp(desc, '\*\.([a-zA-Z0-9_]+)');
            if ~isempty(exm{3})
                newx = desc(exm{3}{1}(1, 1):exm{3}{1}(1, 2));
                xrf_ext(tc + 1, :) = {['*.' newx], desc};
            end
        end
    end
    xrf_ext.x(end, :) = {'*.*', 'All files (*.*)'};
end

% default: no files
fnames = {};

% get or save
if nargin < 5 || isempty(saveas) || (~islogical(saveas) && ~ischar(saveas))
    saveas = false;
end

% request type and/or multiple files
if islogical(saveas) && saveas(1)
    mf = 'off';
    if ischar(nrf)
        dtitle = nrf(:)';
    else
        dtitle = 'Please select or type a filename for saving...';
    end
elseif ischar(saveas)
    mf = 'off';
    dtitle = saveas(:)';
    saveas = false;
elseif nrf > 1
    mf = 'on';
    dtitle = 'Please choose at least one file...';
else
    mf = 'off';
    dtitle = 'Please choose a file...';
end

% parse pattern
[pp, pf, pe] = fileparts(pattern);
if ~isempty(pp)
    op = pwd;
    try
        cd(pp);
    catch xfferror
        neuroelf_lasterr(xfferror);
        pp = '';
    end
end
if ~isempty(pe)
    pe(1) = [];
end
if isempty(pe)
    pe = '*';
end

% build pattern list
pattlist = cell(0, 2);
useall = false;

% get valid extensions
exl = fieldnames(ext);
if numel(pe) > 1 && pe(1) == '(' && pe(end) == ')'
    pe = xrf_ext.n.splittocellc(pe(2:end-1), '|');
    pattlist = cell(numel(pe), 2);
    for fc = 1:numel(pe)
        exf = find(strcmpi(exl, pe{fc}));
        if ~isempty(exf)
            fmt = ext.(exl{exf(1)});
            spec = formats.(lower(fmt{1}(end-2:end)))(fmt{2});
            exf = find(strcmpi(spec.Extensions, pe{fc}));
            desc = spec.Description{exf(1)};
            if isempty(strfind(desc, '*.'))
                pattlist(fc, :) = {['*.' pe{fc}], [desc ' (*.' pe{fc} ')']};
            else
                [exm{1:3}] = regexp(desc, '\*\.([a-zA-Z0-9_]+)');
                if ~isempty(exm{3})
                    newx = desc(exm{3}{1}(1, 1):exm{3}{1}(1, 2));
                    pattlist(fc, :) = {['*.' newx], desc};
                end
            end
        end
    end
else
    exf = find(strcmpi(exl, pe));
    if ~isempty(exf)
        fmt = ext.(exl{exf(1)});
        spec = formats.(lower(fmt{1}(end-2:end)))(fmt{2});
        exf = find(strcmpi(spec.Extensions, pe));
        desc = spec.Description{exf(1)};
        if isempty(strfind(desc, '*.'))
            pattlist(end + 1, :) = {['*.' pe], [desc ' (*.' pe ')']};
        else
            [exm{1:3}] = regexp(desc, '\*\.([a-zA-Z0-9_]+)');
            if ~isempty(exm{3})
                newx = desc(exm{3}{1}(1, 1):exm{3}{1}(1, 2));
                pattlist(end + 1, :) = {['*.' newx], desc};
            end
        end
    end
end

% add all types to pattern list for defined types
if ~strcmp(pe, '*')
    pattlist(end + 1, :) = {'*.*', 'All files (*.*)'};

% otherwise create list !
else
    pattlist = [xrf_ext.x(xrf_ext.l, :); xrf_ext.x(~xrf_ext.l, :)];
    useall = true;
end

% check for filename
if strcmp(pf, '*')
    pf = cell(0, 1);
else
    pf = {[pf '.' pe]};
end

% try getting file
try
    if ~saveas
        if xrf_ext.n.mainver() > 6
            [fnames, fpath, fidx] = uigetfile(pattlist, dtitle, pf{:}, 'MultiSelect', mf);
            if useall && numel(fidx) == 1 && fidx > 0
                fidi = [find(xrf_ext.l), find(~xrf_ext.l)];
                xrf_ext.l(:) = false;
                xrf_ext.l(fidi(fidx)) = true;
            end
        else
            [fnames, fpath] = uigetfile(pattlist, dtitle, pf{:});
        end
    else
        [fnames, fpath] = uiputfile(pattlist, dtitle, pf{:});
    end
catch xfferror
    oeo = xfferror;
    if ~isempty(pp)
        try
            cd(op);
        catch xfferror
            neuroelf_lasterr(xfferror);
        end
    end
    warning('neuroelf:xff:callError', ...
        'Error calling uigetfile/uiputfile(...): %s.', oeo.message);
    return;
end
if ~isempty(pp)
    try
        cd(op);
    catch xfferror
        neuroelf_lasterr(xfferror);
    end
end

% return on cancel
if isequal(fnames, 0) || isequal(fpath, 0)
    fnames = {};
    return;
end

% make cell array
if ~iscell(fnames)
    fnames = {fnames};
end

% prepend path
fpath = strrep(fpath, filesep, '/');
if ~isempty(fpath) && fpath(end) ~= '/'
    fpath(end + 1) = '/';
end
for fc = 1:numel(fnames)
    fnames{fc} = [fpath fnames{fc}];
end
