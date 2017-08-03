function outmat = normsound(inmat);

inmat = inmat + abs(min(inmat)); %makes it definitely above zero
therange = range(inmat);
outmat = inmat ./therange;
outmat = outmat .*2 -1;
return

