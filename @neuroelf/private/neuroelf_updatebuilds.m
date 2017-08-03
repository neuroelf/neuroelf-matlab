function neuroelf_updatebuilds
% neuroelf_updatebuilds  - update build info in file headers
%
% FORMAT:       neuroelf_updatebuilds
%
% No input/output fields.
%
% Note: this function should only be used on Linux/Mac machines with the
%       "touch" command line utility
%
%       also, please be aware that the "Date: " field will NOT be udpated!

% Version:  v1.0
% Build:    14091216
% Date:     Sep-12 2014, 4:13 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2013, 2014, Jochen Weber
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

% no arguments, no check, get build information on files
dispp('NeuroElf code maintenance -> updating build information in files...');
dispp(' -> reading build information (from *.m/c/cpp/h/?ff/tfg files)...');
[build, builds, bvers] = neuroelf_build;

% remove entries that are not to be updated
builds(~cellfun('isempty', regexpi(builds(:, 1), '^_files')), :) = [];
builds(cellfun(@isnan, builds(:, 2)), :) = [];

% test whether build concurs with date in filesystem
nep = [neuroelf_path, '/'];
for fc = 1:size(builds, 1)

    % get dir information on file
    nefname = [nep, builds{fc, 1}];
    dinfo = dir(nefname);

    % parse into build number
    dbuild = str2double(datestr(dinfo.datenum, 'yymmddHH'));
    ddate = strrep(datestr(dinfo.datenum, 'mmm-dd yyyy, HH:MM PM'), '  ', ' ');
    dyear = datestr(dinfo.datenum, 'yyyy');

    % must be within an hour
    if abs(builds{fc, 2} - dbuild) > 1
        dispp(sprintf(' -> updating %s from Build: %d to %d.', ...
            builds{fc, 1}, builds{fc, 2}, dbuild));

        % read file
        nefcont = asciiread(nefname);
        newcont = nefcont;

        % update version
        if isempty(strfind(nefcont, ['Version:  v' bvers]))
            newcont = regexprep(newcont, 'Version:  v\d+\.\d[a-z]?', ...
                ['Version:  v' bvers]);
        end

        % replace Build
        newcont = regexprep(newcont, ...
            sprintf('Build:(\\s+)%d', builds{fc, 2}), ...
            sprintf('Build:$1%d', dbuild));

        % compile date
        newcont = regexprep(newcont, ...
            'Date:\s+[A-Z][a-z][a-z]\-\d\d?\s+20\d\d,\s+\d\d?:\d\d\s?[aApP][mM](\sEST)+', ...
            ['Date:     ' ddate ' EST']);

        % add year
        if isempty(regexpi(newcont, ...
            ['copyright\s\(c\)\s\d+[ ,\-0-9]+' dyear ',\s+jochen\sweber'])) && ...
            isempty(regexpi(newcont, ...
            ['copyright\s\(c\)\s' dyear ',\s+jochen\sweber']))
            newcont = regexprep(newcont, ...
                '(Copyright\s+\(c\)\s\d+[ ,\-0-9]+)(Jochen\sWeber)', ...
                ['$1' dyear ', $2']);
        end

        % not updated?
        if strcmp(nefcont, newcont)
            dispp('   -> failed.');
            continue;
        end

        % write back
        asciiwrite(nefname, newcont);

        % and touch
        tbuild = datestr(dinfo.datenum, 'yymmddHHMM.SS');
        try
            [ss, touched] = system(sprintf('touch -t %s %s', tbuild, nefname));
            if ss ~= 0
                dispp(['   -> failed to touch ' nefname ': ' touched]);
            end
        catch ne_eo;
            dispp(['   -> failed to touch ' nefname ': ' ne_eo.message]);
        end
    end
end
