function wrpembed
%WRPEMBED Embed (or disembed) a binary file in a picture (BMP or PNG).
%   WRPEMBED will ask for an input file, and a carrier file. If the carrier
%   file is not specified, WRPEMBED will attempt to disembed a stream from
%   the carrier (first input). Next, WRPEMBED will request a numeric PIN to
%   be used to seed a random number generator stream (RandStream, using the
%   MT19937AR algorithm), and then a pass phrase, which is used to further
%   obfuscate the binary stream.

% Version:  v1.1
% Build:    16122114
% Date:     Dec-21 2016, 2:30 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

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

% request first file
[i1, ip] = uigetfile({'*.png'; '*.bmp'; '*.*'}, 'Please select file to embed/disembed from');
if isequal(i1, 0)
    return;
elseif ~isempty(ip)
    cd(ip);
end

% request second file
[n2, ip] = uigetfile({'*.png'; '*.bmp'}, 'Please select carrier (for embedding, or cancel)');
try
    if isequal(n2, 0)
        i1 = imread(i1);
        i1 = i1(:);
        l1 = floor(numel(i1) / 8);
        i2 = [];
    else
        n1 = i1;
        nn = numel(n1);
        f1 = fopen(n1, 'r');
        i1 = fread(f1, [1, Inf])';
        l1 = numel(i1);
        fclose(f1);
        if ~isempty(ip)
            cd(ip);
        end
        i2 = imread(n2);
        s2 = size(i2);
        l2 = floor(numel(i2) / 8);
    end
catch ie
    rethrow(ie);
end

% request PIN and generate random number stream
ip = inputdlg('PIN:', 'Input', 1, {'0'});
if isempty(ip) || isempty(ip{1})
    return
end
try
    rs = RandStream('mt19937ar', 'Seed', ceil(str2double(ip{1})));
catch ie
    rethrow(ie);
end

% request pass phrase
ip = inputdlg('Pass phrase:', 'Input', 1, {''});
if isempty(ip)
    return;
else
    ip = ip{1};
    if isempty(ip)
        ip = '  ';
    end
    ip = mod(cumsum(89 .* double(ip(:))), 256);
end

% disembed
if isempty(i2)

    % first get the least significant bit of all pixels
    i1 = sum(reshape(repmat(2.^((0:7)'), l1, 1) .* mod(double(i1(1:8*l1)), 2), 8, l1), 1)';
    i1 = wrpe(i1, ip, @minus, rs);

    % extract filename length, filename, file length and content
    nn = 256 * i1(1) + i1(2);
    if nn > (numel(i1) - 2)
        error('Error disembedding');
    end
    n1 = i1(3:2+nn)';
    if any(n1 < 32 | n1 > 127)
        error('Error disembedding');
    end
    l1 = sum(i1(3+nn:6+nn) .* (2 .^ (24:-8:0)'));
    if l1 > (numel(i1) - (6 + nn))
        error('Error disembedding');
    end
    i1 = i1(7+nn:6+nn+l1);

    % save file
    [op, on, oe] = fileparts(char(n1));
    [n1, ip] = uiputfile(['*', oe], 'Save outputfile...', char([on, oe]));
    if isequal(n1, 0)
        return;
    end
    if ~isempty(ip)
        cd(ip);
    end
    f1 = fopen(n1, 'w');
    fwrite(f1, i1(:), 'uint8');
    fclose(f1);

% embed
else

    % pack filename length, filename, file length and content
    nn = [floor(nn / 256), mod(nn, 256)];
    ln = mod(floor(l1 ./ (2 .^ (24:-8:0))), 256);
    i1 = [nn, double(n1), ln, i1(:)'];

    % add 0s
    if numel(i1) < l2
        i1(l2) = 0;

    % resize image
    elseif numel(i1) > l2
        rf = sqrt(numel(i1) / l2);
        i2 = imresize(i2, ceil(rf .* s2(1:2)), 'cubic');
        n2 = ['r', n2];
        s2 = size(i2);
        l2 = floor(numel(i2) / 8);
        if numel(i1) < l2
            i1(l2) = 0;
        end
    end

    % encrypt
    i1 = wrpe(i1, ip, @plus, rs);

    % clear least significant bit
    i2 = bitand(i2(:), uint8(254));

    % embed in image
    i1 = uint8(bitand(uint8(repmat(i1(:)', 8, 1)), repmat(uint8(2 .^ (0:7)'), 1, l2)) > 0);

    % combine and save
    i2(1:numel(i1)) = i2(1:numel(i1)) + i1(:);
    imwrite(reshape(i2, s2), n2);
end

% sub-function for de/enrcyption
function m = wrpe(m, k, f, rs)

% expand raw pass phrase
nm = numel(m);
k = repmat(k(:), ceil(nm / numel(k)), 1);

% process message
m = m(:);
m = mod(f(f(double(m), k(1:nm)), ceil(256 .* rs.rand(nm, 1))), 256);
