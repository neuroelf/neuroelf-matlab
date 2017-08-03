function y = str2(x,newsize);

if nargin == 1, newsize ==2; end;
	
% returns a digit string of size newsize from a
% number of size x, padding with 0s when necessary
y = num2str(x);
s1 = length(y);
for i = 1: (newsize-s1)
  y = ['0' y];
end
