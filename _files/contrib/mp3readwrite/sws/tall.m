function x = tall(y)
if length(size(y)) > 2 
  error('Too many dims for tall.')
end

[nr nc] = size(y);
x = nr;
