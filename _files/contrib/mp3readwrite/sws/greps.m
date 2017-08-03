function y = greps(datum,inchar,outchar);

datum(find(datum == inchar)) = outchar;
y = datum;




