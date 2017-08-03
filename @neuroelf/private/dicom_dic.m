function dic = dicom_dic(tag, dchoice)
% dicom_dic  - resolve DICOM tag
%
% FORMAT:       tag = dicom_dic(key)
%   - OR -      key = dicom_dic(tag)
%   - OR -      dic = dicom_dic
%
% Input fields:
%
%       key         if given, must be a 'xxxx.xxxx' key
%       tag         if given, must be a valid tag name
%
% Output fields:
%
%       key, tag    see above
%       dic         full dictionary structure
%
% Note: the main dictionary itself is taken from the DCMTK
%       application suite provided by OFFIS e.v. available
%       at http://dicom.offis.de/

% Version:  v0.9b
% Build:    11050712
% Date:     Apr-09 2011, 2:11 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2011, Jochen Weber
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

% persistent variable
persistent dcm_dic;
if isempty(dcm_dic)

    % set main dictionary
    dcm_dic = struct;

    % find dictionaries
    dictfiles = findfiles([neuroelf_path('dicom') filesep '*' filesep], '*.dic', ...
        'oneperdir=1');

    % parse each dictionary
    for dc = 1:length(dictfiles)

        % get dictionary
        try
            ddict = dicom_parsedict(dictfiles{dc});
        catch ne_eo;
            error( ...
                'xff:BadFileContent', ...
                'Error parsing DICOM dictionary: %s (%s).', ...
                dictfiles{dc}, ne_eo.message ...
            );
        end

        % put into array
        dname = fileparts(dictfiles{dc});
        [dnamec{1:2}] = fileparts(dname);
        dcm_dic.(dnamec{2}) = ddict;
    end
end

% dict choice
if nargin > 1 && ...
    ischar(dchoice) && ...
   ~isempty(dchoice) && ...
    isfield(dcm_dic, dchoice(:)')
    dchoice = dcm_dic.(dchoice(:)');
else
    dchoice = dcm_dic.OFFIS;
end

% no valid tag, return complete dict
if nargin < 1 || ...
   ~ischar(tag) || ...
    isempty(tag)

    % special cases
    if nargin > 0 && ...
        isa(tag, 'double') && ...
        numel(tag) == 1 && ...
        round(tag) > 0 && ...
        tag < 6

        switch (round(tag))
            case {1}
                dic = dchoice.Dictionary;
            case {2}
                dic = dchoice.KeyToTag;
            case {3}
                dic = dchoice.TagToKey;
            case {4}
                dic = dcm_dic.OFFIS;
            otherwise
                dic = fieldnames(dcm_dic);
        end
        return;
    end
    dic = dchoice.TagToKey;

% key
elseif ~isempty(regexpi(tag(:)', '^[0-9a-f]{4}[^0-9a-z]?[0-9a-f]{4}$'))

    % get key
    tag = tag(:)';
    if length(tag) == 9
        tag = [tag(1:4) tag(6:9)];
    end
    tag = ['k_' upper(tag(1:4)) '_' upper(tag(5:8))];

    % check if key exists in dict
    if isfield(dchoice.Dictionary, tag)

        % return entry
        dic = dchoice.Dictionary.(tag);

    % error otherwise
    else
        error( ...
            'xff:BadArgument', ...
            'Unknown DICOM key: %s.', ...
            [tag(3:6) '.' tag(8:11)] ...
        );
    end

% label key
elseif strcmp(makelabel(tag(:)'), tag(:)')

    % check if tag exists
    tag = tag(:)';
    if isfield(dchoice.TagToKey, tag)

        % return entry
        dic = dchoice.Dictionary.(dchoice.TagToKey.(tag));

    % error otherwise
    else
        error( ...
            'xff:BadArgument', ...
            'Unknown DICOM tag: %s.', ...
            tag ...
        );
    end

% otherwise
else
    error( ...
        'xff:BadArgument', ...
        'Invalid tag/key format given.' ...
    );
end


% %%%  subfunctions


function ddic = dicom_parsedict(dictfile)

    % read dictionary
    rawdict = splittocell(asciiread(dictfile), char([10, 13]), 1, 1);

    % init dictionary
    ddic = struct;
    dct = struct;
    ktt = struct;
    ttk = struct;

    % parse dictionary lines
    for lc = 1:length(rawdict)

        % reject comments
        if isempty(rawdict{lc}) || ...
            rawdict{lc}(1) == '#'
            continue;
        end

        % and malformatted lines
        line = splittocell(rawdict{lc}, char([9, 32]), 1, 1);
        if length(line) ~= 5 || ...
            length(line{1}) ~= 11 || ...
            line{1}(1) ~= '(' || ...
            line{1}(end) ~= ')'
            continue;
        end

        % extract key and name
        dkey = ['k_' upper(line{1}(2:5)) '_' upper(line{1}(7:10))];
        dtag = line{3};

        % check key and name
        if ~strcmp(dkey, makelabel(dkey)) || ...
           ~strcmp(dtag, makelabel(dtag))
            continue;
        end

        % create entries
        dstr = struct( ...
            'key', dkey, ...
            'vr',  line{2}, ...
            'tag', dtag, ...
            'vm',  line{4}, ...
            'ver', line{5});

        % set entries only once!
        if ~isfield(dct, dkey)
            dct.(dkey) = dstr;
            ktt.(dkey) = dtag;
            ttk.(dtag) = dkey;
        end
    end

    % put structs in main struct
    ddic.Dictionary = dct;
    ddic.KeyToTag = ktt;
    ddic.TagToKey = ttk;

% end of function ddic = parsedict(dictfile)
