function s = slydetrend(s)
% slydetrend  - slyly detrends signal

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

ss = size(s);
if numel(s) == length(s)
    rs = [1,1];
    s = s(:);
    sn = numel(s);
else
    rs = [1, prod(ss(2:end))];
    s = reshape(s, [ss(1), rs(2)]);
    sn = ss(1);
end
sd = sn - 1;
fs = reshape(s(1,:), rs);
ls = reshape(s(end,:), rs);

ps = [fs;fs;s;ls;ls];
ps = convn(ps, [1;1;1]./3, 'same');
ps = diff(ps(2:end-2,:));
mdp = mean(ps);
sdp = std(ps);
gdi = abs((ps-ones(sn, 1)*mdp))<2*(ones(sn, 1)*sdp);
ps(~gdi) = 0;
mdp = sum(ps) ./ sum(gdi);
s = s - ([0:sd]' * mdp) + ((sd/2).*(ones(sn,1)*mdp));
reshape(s, ss);
