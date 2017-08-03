function y = fileext(x);

len = length(x);
if len < 1, error('No data.'); end

thischar = len;
y = [];

while (thischar > 1) & (x(1,thischar) ~= '.')
y = [x(1,thischar) y];
thischar = thischar -1;
end
if x(1,thischar) ~= '.',  y = []; end;

return
