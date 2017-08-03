function samps=msecs2samps(msecs,srate)

% sample <-> time conversion
samps=floor(msecs*srate/1000);
