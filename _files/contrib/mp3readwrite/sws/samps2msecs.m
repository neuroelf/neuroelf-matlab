function msecs=samps2msecs(samps,srate)

% sample <-> time conversion
if isequal(samps,0)
	msecs=1000*samps/srate;
else
	msecs = 1000*(samps)/srate;
end
