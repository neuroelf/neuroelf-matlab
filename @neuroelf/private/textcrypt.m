function t = textcrypt(t, k, d)
% textcrypt  - en/decrypt a text message with a rotating key matrix
%
% FORMAT:       t = textcrypt(t, key [, decrypt])
%
% Input fields:
%
%       t           text to en- or decrypt
%       key         key phrase (up to 48 significant characters)
%       decrypt     decrypt flag (default: false)
%
% Output fields:
%
%       t           encrypted (or decrypted) text

% Version:  v0.9a
% Build:    10051716
% Date:     May-17 2010, 10:48 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, Jochen Weber
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

% persistent initial key
persistent tcki;
if isempty(tcki)
    tcki = double([ ...
        '"i@o><K8e''/-HyLX$URE4c|3b2j}x*v&wd%6`T^7q!thn.N#aWg]?mIz)=fk;9uplFQr1O([sYA~0SJM\_+G5B{PDC,ZV:'; ...
        'keK<h#La97d^$5xyf|VAI[2%B8;EjTX(H*bP3/Sc04N`pzg,)!n~l1v@oQqsGwJ{_''u\?.YtFW}=-ZrR6+]:DCiU>m"M&O'; ...
        '-SKw6@yIAb>#D{Qme79GsVUOc:(F&\XEapCZ0<)]Win''hHT.fBJdM$|^j"*?%Lo`P8;uNv,2t~/r1R5}Y[=gz_+k43lxq!'; ...
        '/X?%hg[sR<lK\up_"&|$yjtVA3>b`@)qO{MQ=F#G9Zc!*~Lf5,]^;z-EH10C}oiTr+m(kUN''n2xBdJW6D4e7vYPI:8Saw.'; ...
        'q@\W{wOCx.J%b5~!92zghG;6j]S"3=&`fI''LE/U}KsA|<dNM#DcnR0tV-H^T17iaom,YF_B4ylv:$Z?u+e(P>rQp*8[kX)'; ...
        ':{hiTG#Wr9myMuN!s)<?(OgKVp4kB|DS\l1jd&[Lt]^3-6''aAo_nxU"f*~X;5YP.bz8%eqZ+=@J7Ev,$RIwC>F`Q/c20H}'; ...
        '?OA,p^-`sXNz#C2.&{ME+]LvlDk|mZ%JU@86_/\*tYW)i=F[BwaPS7Toj}0KRu1x$e4gI>(;drQ"''39nbGyf:h~5HqV!c<'; ...
        '47W\[neDx;dB5MV{&qR1y-cJi]6HX`*$%ATwm^P~a<K+.u9}Ut/FEg''0rYl@N,:p>"Lz?#Zb3jf(ov)2|I8h=O!kGsC_QS'])' - 32;
end

% argument check
if nargin < 2 || ...
   ~ischar(t) || ...
    numel(t) ~= size(t, 2) || ...
   ~ischar(k) || ...
    isempty(k)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin < 3 || ...
   ~islogical(d) || ...
    numel(d) ~= 1
    d = false;
end

% create numeric key (with minimum length)
k = mod(double(k(:)') - 33, 94) + 1;
kl = numel(k);
if kl < 16
    k = repmat(k, 1, ceil(16 / kl));
    k = k(1:16);
    kl = 16;
end
while kl > 48
    if mod(kl, 2) ~= 0
        k(end + 1) = k(1);
    end
    k = mod(k(1:2:end) + k(1:2:end) .* k(2:2:end) - k(end:-2:2), 94) + 1;
end

% initial key indices
ki = tcki;

% which initial key number and snippet size
kp = true(1, 94);
kn = mod(sum(k), 8) + 1;
ss = mod(sum(k), 32) + 63;

% shuffle keys with key
for sc = 1:8
    kp(:) = true;
    k = mod(k + ki(k + 94 * mod(k * 5 + (1:kl), 8)), 94) + 1;
    kp(k) = false;
    ki = [ki(~kp, :); ki(kp, :)];
    ks = sum(k .* k);
    [ko, koi] = sort(ki(mod(ks:ks+7, 87) + 1, mod(ks, 8) + 1));
    ki = ki(:, koi);
end
ki = ki';

% add marker for OK decryption
if ~d
    t(end+1:end+10) = '$TEXTCRYPT';
end

% find regular characters
tr = (t > 32 & t < 127);

% get raw values
r = double(t(tr)) - 32;

% initialize counters
cp = 1;
nr = numel(r);

% decrypt
if d

    % until everything is decrypted
    while cp <= nr

        % get encrypted snippet
        ss = min(ss, nr - cp + 1);
        rs = r(cp:cp + ss - 1);

        % and get corrected version
        rc = mod(rs - (2 + ki(kn, 1:ss)), 94) + 1;

        % and determine next key and snippet size
        knn = mod(sum(rs), 8) + 1;
        ssn = mod(sum(rs), 32) + 63;

        % get decoding matrix
        [dv, dc] = sort(ki(kn, :));

        % then patch key first
        kr = mod(rs(end) + 13 * [rs(1):94, 1:rs(1)-1], 94) + 1;
        ki(kn, :) = ki(kn, kr);

        % and then decode snippet
        rs = dc(rc);

        % and shuffle new key
        ki(knn, :) = ki(knn, ki(kn, :));
        kn = knn;

        % go on
        r(cp:cp + ss - 1) = rs;
        cp = cp + ss;
        ss = ssn;
    end

% encrypt
else

    % until everything is crypted
    while cp <= nr

        % get next snippet
        ss = min(ss, nr - cp + 1);
        rs = r(cp:cp + ss - 1);

        % get encoded version
        rs = mod(ki(kn, rs) + ki(kn, 1:ss), 94) + 1;

        % then patch key
        kr = mod(rs(end) + 13 * [rs(1):94, 1:rs(1)-1], 94) + 1;
        ki(kn, :) = ki(kn, kr);

        % go on
        r(cp:cp + ss - 1) = rs;
        cp = cp + ss;

        % and determine next key and snippet size
        knn = mod(sum(rs), 8) + 1;
        ss = mod(sum(rs), 32) + 63;

        % and shuffle new key
        ki(knn, :) = ki(knn, ki(kn, :));
        kn = knn;
    end
end

% then retransfer to char
t(tr) = char(r + 32);

% marker found
if d
    if ~strcmp(t(end-9:end), '$TEXTCRYPT')
        error( ...
            'neuroelf:InvalidPassword', ...
            'Incorrect password provided.' ...
        );
    end
    t(end-9:end) = [];
end
