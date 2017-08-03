function [string, number] = parseline(x);

theline = x;
%theline = 'abcde 0.1234';

L = length(theline);

positn = 1;
str1 = [];

%Go through the line character-by-character until white space

while [1]
  tempositn = theline(1, positn);
  if isspace(tempositn) 
    positn = positn+ 1;
	break; end;
  str1 = [str1 tempositn];
  positn = positn + 1;
end 

string = str1;

%Continue after white space until end of line

str1 = [];
startc = positn;
for positn = startc:L
  tempositn = theline(1, positn);
  str1 = [str1 tempositn];
end 

number = str2num(str1);



