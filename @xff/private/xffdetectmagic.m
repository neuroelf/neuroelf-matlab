function detmag = xffdetectmagic(filename, mag)
% xff::_detectmagic  - detect filetype from magic token
%
% THIS FUNCTION IS AN INTERNAL FUNCTION OF THE CLASS
%
% @xff
%
% AND IT SHOULD NEVER BE CALLED FROM OUTSIDE THE CLASS CODE

% Version:  v1.0
% Build:    16012916
% Date:     Jan-29 2016, 4:09 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2016, Jochen Weber
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

% default: no magic
detmag = '';

% read file first
try
    fileinfo = dir(filename);
    if numel(fileinfo) ~= 1 || fileinfo.isdir
        return;
    end
    filesize = fileinfo.bytes;
    filecont = transio(filename, 'ieee-le', 'uint8', 0, [1, filesize]);
    firstblk = filecont(1, 1:min(2048, filesize));
catch xfferror
    neuroelf_lasterr(xfferror);
    return;
end

% try magics
for ftc = 1:numel(mag)

    % get name, range, type, token, and infer detected type
    mmag = mag(ftc);
    mname  = lower(mmag.name);
    mnames = find(mname == '_');
    if isempty(mnames) || mnames(1) < 4
        warning('neuroelf:xff:invalidToken', 'Invalid Magic token tag: %s.', mname);
        continue;
    end
    mrange = mmag.range;
    mtype  = lower(mmag.type);
    mtoken = mmag.magic;
    mdtype = mname(1:mnames(1)-1);
    if mrange(2) > filesize
        if strcmp(mtype, 'hex')
            continue;
        else
            mrange(2) = filesize;
        end
    end
    if mrange(2) <= 2048
        mdcont = char(firstblk(1, mrange(1):mrange(2)));
    else
        mdcont = char(filecont(1, mrange(1):mrange(2)));
    end

    % default: no match
    mmatched = false;

    % what type
    switch (mtype)
        case 'hex'
            if numel(mtoken) == numel(mdcont) && all(double(mtoken(:)') == double(mdcont(:)'))
                mmatched = true;
            end
        case 'regexp'
            if ~isempty(regexp(char(mdcont), mtoken, 'once'))
                mmatched = true;
            end
        case 'regexpi'
            if ~isempty(regexpi(char(mdcont), mtoken))
                mmatched = true;
            end
        case 'strfind'
            if ~isempty(strfind(char(mdcont), mtoken))
                mmatched = true;
            end
    end

    % leave loop ?
    if mmatched
        detmag = mdtype;
        break;
    end
end
