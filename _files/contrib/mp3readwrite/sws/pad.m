function y = pad(x,n);	%this function adds extra spaces to text lines
						%to create matrices with equal dimensions
l = length(x);

spacesneeded = n-l;

for i = 1:spacesneeded
x = [x ' '];
end

y = x;


