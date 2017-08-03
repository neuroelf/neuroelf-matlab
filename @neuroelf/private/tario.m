function [varargout] = tario(varargin)
% tario  - read TAR files
%
% FORMAT:       tarobject = tario(filename)
%
% Input fields:
%
%       filename    filename of TAR file to read
%
% Output fields:
%
%       tarobject   xff object
%
% See also xff.

% See http://www.gnu.org/software/tar/manual/html_node/Standard.html

% Version:  v1.1
% Build:    16060314
% Date:     Jun-03 2016, 2:30 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2016, Jochen Weber
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

% persistent VR dict and empty file
persistent my_newtar;
if isempty(my_newtar)
    newtar = xff('new:tar');
    my_newtar = getcont(newtar);
    delete(newtar);
end

% argument check
if nargin < 1 || ~ischar(varargin{1}) || isempty(varargin{1})
    error('neuroelf:general:badArgument', 'Bad or missing argument for tario.');
end
filename = varargin{1}(:)';
if exist(filename(:)', 'file') ~= 2
    error('neuroelf:general:badArgument', 'TAR file doesn''t exist.');
end
tarcont = my_newtar;

% try to open file
try
    fid = fopen(filename, 'r');
    if fid < 1
        error('neuroelf:fileIO:fileOpenError', 'Cannot open TAR file for reading.');
    end
catch ne_eo;
    rethrow(ne_eo);
end
fseek(fid, 0, 1);
fsize = ftell(fid);
fseek(fid, 0, -1);

% start with a modest list size
names = cell(4096, 1);
fhpos = zeros(4096, 1);
fcpos = zeros(4096, 1);
fclen = zeros(4096, 1);
types = char(48 * ones(4096, 1));

% start reading until EOF reached
fc = 0;
rpos = 1;
mpos = fsize - 510;
ov11 = 8 .^ (10:-1:0);
ov12 = 8 .^ (11:-1:0);
while rpos < mpos

    % increase file counter
    fc = fc + 1;

    % current position
    fhpos(fc) = rpos;

    % read one header segment
    fh = fread(fid, [1, 512], 'uint8=>uint8');
    rpos = rpos + 512;

    % get filetype
    ft = fh(157);
    if ft == 0
        ft = 48;
    end

    % filecontent length
    if fh(136) > 32
        fcl = sum(ov12 .* double(fh(125:136) - 48));
    else
        fcl = sum(ov11 .* double(fh(125:135) - 48));
    end

    % extended header type 'L' (@LongLink) or 'x' (PaxHeader/)
    fnh = [];
    if any(ft == 'Lx') && fcl > 0

        % re-read name
        fnrlen = 512 * ceil(fcl / 512);
        fnh = fread(fid, [1, fnrlen], 'uint8=>uint8');
        rpos = rpos + fnrlen;
        fnh = fnh(1:fcl);
        if fnh(end) == 0
            fnh(end) = [];
        end
        if ft == 'x'
            fnbeg = findfirst(fnh == 32);
            fnend = str2double(char(fnh(1:fnbeg-1)));
            fnh = fnh(fnbeg+1:fnend);
            if numel(fnh) > 5 && all(fnh(1:5) == [112, 97, 116, 104, 61])
                fnh(1:5) = [];
            end
            if any(fnh(end) == [0, 10, 32])
                fnh(end) = [];
            end
        end
        fhpos(fc) = rpos;
        fh = fread(fid, [1, 512], 'uint8=>uint8');
        rpos = rpos + fnrlen;
        ft = fh(157);
        if ft == 0
            ft = 48;
        end
        if fh(136) > 32
            fcl = sum(ov12 .* double(fh(125:136) - 48));
        else
            fcl = sum(ov11 .* double(fh(125:135) - 48));
        end
    end

    % content
    fcpos(fc) = rpos;

    % filename length
    fnl = findfirst(fh == 0);

    % end of archive
    if fnl == 1 && all(fh == 0)
        fc = fc - 1;
        break;
    end

    % store name
    if isempty(fnh)
        if fnl > 100
            fnl = 101;
        end
        names{fc} = char(fh(1:fnl-1));
    else
        names{fc} = char(fnh);
    end

    % store length
    fclen(fc) = fcl;

    % increase position
    if mod(fcl, 512) ~= 0
        fcl = 512 * ceil(fcl / 512);
    end
    fseek(fid, fcl, 0);
    rpos = rpos + fcl;
    
    % store type
    types(fc) = ft;

    % extend lists
    if mod(fc, 4096) == 0
        names(end+4096) = {[]};
        fhpos(end+4096) = 0;
        fcpos(end+4096) = 0;
        fclen(end+4096) = 0;
        types(end+1:end+4096) = '0';
    end
end

% close file
fclose(fid);

% store data
if numel(names) > fc
    names(fc+1:end) = [];
    fhpos(fc+1:end) = [];
    fcpos(fc+1:end) = [];
    fclen(fc+1:end) = [];
    types(fc+1:end) = [];
end
tarcont.Name = names;
tarcont.Type = types;
tarcont.HeadPos = fhpos;
tarcont.ContPos = fcpos;
tarcont.ContLen = fclen;

% set
hfile = xff('new:tar');
tarcont.RunTimeVars.xffID = hfile.RunTimeVars.xffID;
setcont(hfile, tarcont);

% return
varargout{1} = hfile;
