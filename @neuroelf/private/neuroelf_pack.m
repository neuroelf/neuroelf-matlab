function neuroelf_pack(tfolder, update)
% neuroelf_pack  - create installation/update/diff package
%
% FORMAT:       neuroelf_pack(tfolder [, update])
%
% Input fields:
%
%       tfolder     target folder name
%       update      optional string, if given only files newer than date
%
% No output fields.
%
% Note: this function uses mpackage to create
%
%  - NeuroElf_vXXfullYY.m and .mat
%  - NeuroElf_vXXYY.m and .mat
%
% Both types of installation contain the necessary files!

% Version:  v1.1
% Build:    17080412
% Date:     Aug-04 2017, 12:45 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2016, 2017, Jochen Weber
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
if nargin < 1 || ~ischar(tfolder) || exist(tfolder(:)', 'dir') ~= 7
    error('neuroelf:general:badArgument', ...
        'First argument must be existing folder name.');
end
if tfolder(end) == filesep
    tfolder(end) = [];
end
if nargin < 2 || ~ischar(update)
    update = [];
else
    if strcmpi(update(:)', 'sinceinstall')
        upds = true;
        update = 'neuroelf';
    else
        upds = false;
    end
    try
        if exist(update(:)', 'file') == 2
            udir = dir(which(update(:)'));
            update = datenum(udir.date) - 1 / 1440;
        else
            update = datenum(update(:)');
        end
        if upds
            update = update + 1 / 96;
        end
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        update = [];
    end
end

% get current folder name
nelffld = neuroelf_path;
[nelfp{1:2}] = fileparts(nelffld);

% update builds in files
if ~isempty(regexpi(computer, '^maci'))

    % touch neuroelf_build first
    system(sprintf('touch %s/neuroelf_build.m', nelffld));
    rehash;
    clear neuroelf_build;

    % then update build information in files
    neuroelf_updatebuilds;
    neb = neuroelf_build;

    % next re-create the splash text images
    disp(' -> re-creating splash text images...');
    neuroelf_sptexts;

    % and finally clean up a bit prior to making the packages
    disp(' -> temporarily moving binaries out of the way...');
    nelftp = [nelfp{1} filesep 'NE_files'];
    if exist(nelftp, 'dir') ~= 7
        mkadir(nelftp);
    end
    cachefiles = findfiles([nelffld '/_files/cache'], {'*cache*.mat', 'setup*.done'}, '-d1');
    if ~isempty(cachefiles)
        mdelete(cachefiles);
    end
    colinfiles = findfiles([nelffld '/_files/colin'], {'colin_*.*', 'colin.v*'}, '-d1');
    if ~isempty(colinfiles)
        renamefile(colinfiles, nelftp);
    end
    talfiles = findfiles([nelffld '/_files/tal'], 'talairach_ICBMnorm*.*', '-d1');
    if ~isempty(talfiles)
        renamefile(talfiles, nelftp);
    end
    nsynthfiles = findfiles([nelffld '/_files/neurosynth/terms'], '*nii.gz', '-d1');
    if ~isempty(nsynthfiles)
        renamefile(nsynthfiles, nelftp);
    end
    nsynthfiles = findfiles([nelffld '/_files/neurosynth/rawdata'], 'neurosynth_20*.*', '-d1');
    if ~isempty(nsynthfiles)
        renamefile(nsynthfiles, nelftp);
    end
    shenfiles = findfiles([nelffld '/_files/shenparcel'], {'*.rtv'; '*.srf'; '*.voi'}, '-d1');
    if ~isempty(shenfiles)
        renamefile(shenfiles, nelftp);
    end
    if exist([nelffld '/.git'], 'dir') == 7
        movefile([nelffld '/.git'], [nelftp '/.git']);
    end
    if exist([nelffld '/_files/spm/dartel'], 'dir') == 7
        movefile([nelffld '/_files/spm/dartel'], [nelftp '/dartel']);
    end
    if exist([nelffld '/_files/spm/spm8_dartel_template.mat'], 'file') == 2
        movefile([nelffld '/_files/spm/spm8_dartel_template.mat'], nelftp);
    end

    % creating flat ZIP
    fprintf(' -> creating NeuroElf (%d) flat-folder (ZIP) package...\n', neb);
    cd(nelfp{1});
    invsystem(sprintf('zip -9r %s/%s_flat.zip %s', tfolder, nelfp{2}, nelfp{2}));

    % undo-ing temp moves
    disp(' -> moving binary files back in place...');
    colinfiles = findfiles(nelftp, 'colin*.*', '-d1');
    if ~isempty(colinfiles)
        renamefile(colinfiles, [nelffld '/_files/colin/']);
    end
    talfiles = findfiles(nelftp, 'tal*.*', '-d1');
    if ~isempty(talfiles)
        renamefile(talfiles, [nelffld '/_files/tal/']);
    end
    nsynthfiles = findfiles(nelftp, '*.nii.gz', '-d1');
    if ~isempty(nsynthfiles)
        renamefile(nsynthfiles, [nelffld '/_files/neurosynth/terms']);
    end
    nsynthfiles = findfiles(nelftp, 'neurosynth_20*.*', '-d1');
    if ~isempty(nsynthfiles)
        renamefile(nsynthfiles, [nelffld '/_files/neurosynth/rawdata']);
    end
    shenfiles = findfiles(nelftp, 'shen*.*', '-d1');
    if ~isempty(shenfiles)
        renamefile(shenfiles, [nelffld '/_files/shenparcel']);
    end
    if exist([nelftp '/.git'], 'dir') == 7
        movefile([nelftp '/.git'], [nelffld '/.git']);
    end
    if exist([nelftp '/dartel'], 'dir') == 7
        movefile([nelftp '/dartel'], [nelffld '/_files/spm/dartel']);
    end
    if exist([nelftp '/spm8_dartel_template.mat'], 'file') == 2
        movefile([nelftp '/spm8_dartel_template.mat'], [nelffld '/_files/spm/']);
    end
end

% general options
popt = struct;
popt.addpath = nelfp(2);
popt.banner = [ ...
    ' ', char(10), ...
    '  This package contains NeuroElf v' neuroelf_version ' (build ', ...
        sprintf('%d', neuroelf_build), ').', char([10, 10]), ...
    '  Please visit http://neuroelf.net/ for more information.', char([10, 10]), ...
    '  For bug reports and/or feature requests, please contact the', char(10), ...
    '  primary author and maintainer, Jochen Weber at <jw2661@columbia.edu>'];
popt.csyntax = true;
popt.destref = 'mean';
popt.destrup = 2;
popt.dontask = false;
popt.exclude = { ...
    'cache.mat', ...
    'colin_.*\.(hdr|img|rtv|srf|ssm|tsm|v16|vmr)$', ...
    'colin.*(v16|vmr)$', ...
    'neurosynth_20.*', ...
    'shen_parcel.*', ...
    'spm8_dartel_template.mat', ...
    'talairach_ICBMnorm.*', ...
    'Template.*mm_.*', ...
    '.*_term\.nii\.gz$', ...
    '_todo', ...
    '.*\~$', ...
    '.git' ...
    };
popt.finish = {'neuroelf_setup'; 'try, neuroelf_makefiles ask; end'};
if ~isempty(update)
    popt.maxage = 86400 * (now - update);
end
popt.ierrors = true;
popt.postclr = true;
popt.release = 2006;
popt.savev7 = false;
popt.strrep = {[ ...
    '% All rights reserved.', char(10), ...
    '%', char(10), ...
    '% Redistribution and use in source and binary forms, with or without', char(10), ...
    '% modification, are permitted provided that the following conditions are met:', char(10), ...
    '%     * Redistributions of source code must retain the above copyright', char(10), ...
    '%       notice, this list of conditions and the following disclaimer.', char(10), ...
    '%     * Redistributions in binary form must reproduce the above copyright', char(10), ...
    '%       notice, this list of conditions and the following disclaimer in the', char(10), ...
    '%       documentation and/or other materials provided with the distribution.', char(10), ...
    '%     * Neither the name of Columbia University nor the', char(10), ...
    '%       names of its contributors may be used to endorse or promote products', char(10), ...
    '%       derived from this software without specific prior written permission.', char(10), ...
    '%', char(10), ...
    '% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND', char(10), ...
    '% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED', char(10), ...
    '% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE', char(10), ...
    '% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY', char(10), ...
    '% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES', char(10), ...
    '% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;', char(10), ...
    '% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND', char(10), ...
    '% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT', char(10), ...
    '% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS', char(10), ...
    '% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.' char(10)], ...
    ['% BSDLT73e2054aefd869296deef011' char(10)]};
if ~isempty(update)
    upds = regexprep(datestr(now, 'yymmddHH'), '^0*', '');
    popt.addpath = {};
    popt.finish = [{['try, df = which(mfilename); cd(''..''); delete(df); ' ...
        'delete(strrep(df, ''.m'', ''.mat'')); end']}, popt.finish];
    popt.destref = 'neuroelf';
    popt.destrup = 1;
    popt.reqcond = {['strcmp(neuroelf_version, ''', neuroelf_version ''')']};
    popt.require = {'neuroelf_version'};
    popt.update = true;
else
    popt.update = false;
end

% generate package
try
    target = [tfolder '/' nelfp{2}];
    disp(' -> creating NeuroElf installation package...');
    cpfile([nelffld '/ne_eo.m'], [tfolder '/ne_eo.m']);
    mpackage(target, nelffld, popt);
    if popt.update
        pkgcont = load([target '.mat']);
        pkgcont = pkgcont.pkgcont;
        pkgcont.packs{1}{1} = regexprep(pkgcont.packs{1}{1}, '\d+$', '');
        save([target '.mat'], 'pkgcont', '-v6');
        renamefile([target '.m'], [pkgcont.packs{1}{1} '_up_' upds '.m']);
        renamefile([target '.mat'], [pkgcont.packs{1}{1} '_up_' upds '.mat']);
    end
    if ~isempty(regexpi(computer, '^maci'))
        disp(' -> zipping NeuroElf installation package...');
        cd(tfolder);
        invsystem(sprintf('zip -9r %s_install.zip %s.m %s.mat ne_eo.m', nelfp{2}, nelfp{2}, nelfp{2}));
    end
catch ne_eo;
    rethrow(ne_eo);
end
