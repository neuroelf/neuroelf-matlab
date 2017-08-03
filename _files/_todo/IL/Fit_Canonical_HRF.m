function [hrf, fit, e, param, info] = Fit_Canonical_HRF(tc,TR,Runs,T,p, dur)
% function [hrf, fit, e, param, info] = Fit_Canonical_HRF(tc,TR,Runs,T,p, dur)
%
% Fits GLM using canonical hrf (with option of using time and dispersion derivatives)';
%
% INPUTS:
% 
% tc    - time course
% TR    - time resolution
% Runs  - expermental design
% T     - length of estimated HRF
% p     - Model type
% dur - duration (in seconds) of the event
%
% Options: p=0 - only canonical HRF
%          p=1 - canonical + temporal derivative
%          p=2 - canonical + time and dispersion derivative
% 
% OUTPUTS:
%
% hrf   - estimated hemodynamic response function
% fit   - estimated time course
% e     - residual time course
% param - estimated amplitude, height and width
% info  - struct containing design matrices, beta values etc
%
% Created by Martin Lindquist on 10/02/09

% Programmers' notes:
% ----------------------------
% Edited  12/4/09, Tor Wager and Lauren Atlas : Add input duration option

[len,d] = size(Runs);

X = zeros(len,d*p);

[h, dh, dh2] = CanonicalBasisSet(TR, dur);

for i=1:d
    v = conv(Runs(:,i),h);
    X(:,(i-1)*p+1) = v(1:len);

    if (p>1)
        v = conv(Runs(:,i),dh);
        X(:,(i-1)*p+2) = v(1:len);
    end

    if (p>2)
        v = conv(Runs(:,i),dh2);
        X(:,(i-1)*p+3) = v(1:len);
    end
end
    
X = [(zeros(len,1)+1) X];

b = pinv(X)*tc;
e = tc-X*b;
fit = X*b;

b = reshape(b(2:end),p,d)';
bc = zeros(d,1);

for i=1:d,
    if (p == 1)
        bc(i) = b(i,1);
        H = h;
    elseif (p==2)
        bc(i) = sign(b(i,1))*sqrt((b(i,1))^2 + (b(i,2))^2); 
        H = [h dh];
    elseif (p>2)
        bc(i) = sign(b(i,1))*sqrt((b(i,1))^2 + (b(i,2))^2 + (b(i,3))^2);
        H = [h dh dh2];
    end    

end


hrf = H*b';
param = get_parameters2(H*b',T);

info ={};
info.b = b;
info.bc = bc;
info.X = X;
info.H =H;

end

% END MAIN FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Subfunctions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [h, dh, dh2] = CanonicalBasisSet(TR, dur)

len = round(30/TR); % 30 secs worth of images
xBF.dt = TR;  % time bin length, in seconds; Tor: Might be good to chang from TR to SPM's default, TR/16
xBF.length= len;
xBF.name = 'hrf (with time and dispersion derivatives)';
xBF = spm_get_bf(xBF);

v1 = xBF.bf(1:len,1);
v2 = xBF.bf(1:len,2);
v3 = xBF.bf(1:len,3);

% convolve
% input dur is dur in sec; change to units of xBF.dt (elements)
dur = round(dur ./ xBF.dt ); % secs / (secs/sample) = in samples 
v1 = conv(v1, ones(dur, 1)); v1 = v1(1:len);
v2 = conv(v2, ones(dur, 1)); v2 = v2(1:len);
v3 = conv(v3, ones(dur, 1)); v3 = v3(1:len);

% orthogonalize
h = v1;
dh =  v2 - (v2'*v1/norm(v1)^2).*v1;
dh2 =  v3 - (v3'*v1/norm(v1)^2).*v1 - (v3'*dh/norm(dh)^2).*dh;

% normalize amplitude
h = h./max(h);
dh = dh./max(dh);
dh2 = dh2./max(dh2);

end