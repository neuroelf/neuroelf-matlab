function Y = slinterp(X,F)
% Y = slinterp(X,F)  Simple linear-interpolate X by a factor F
%        Y will have ((size(X)-1)*F)+1 points i.e. no extrapolation
% dpwe@icsi.berkeley.edu  fast, narrow version for SWS

% Do it by rows

sx = prod(size(X));

% Ravel X to a row
X = X(1:sx);
X1 = [X(2:sx),0];

XX = zeros(F, sx);

for i=0:(F-1)
  XX((i+1),:) = ((F-i)/F)*X + (i/F)*X1;
end

% Ravel columns of X for output, discard extrapolation at end
Y = XX(1:((sx-1)*F+1));
