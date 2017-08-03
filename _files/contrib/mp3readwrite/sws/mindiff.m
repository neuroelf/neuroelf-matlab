function x=mindiff(y);

thediff=[y(2:end) - y(1:(end-1))];
x=min(thediff);
