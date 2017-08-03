function X = synthtrax(F,M,T,SR,DUR);

% X = synthtrax(F,M,T,SR,DUR)  Synth regular or irreg-sampled sinusoidal
%	representation. F, M and T define data read from an SWI sine-wave 
%	definition file, specifically the frequency and magnitudes of a set of sine
%   oscillators (one per row) along with the time points for each
%   sample.  X is the resulting audio signal, generated at a sampling
%   rate SR (default: 22050 Hz).  If DUR is present and nonzero, the
%   output is exactly this duration in seconds.
%   1998sep18 dpwe@icsi.berkeley.edu
% $Header: $

if (nargin < 3)
% Default to uniform 10ms sampling
 T = 0.010*[0:(size(F,2)-1)];
else
 T=T ./1000; %Scale for Time specified in milliseconds 
end

if (nargin < 4)
  SR = 22050;
end

if (nargin < 5)
  DUR = 0;
end

if DUR == 0
  DUR = max(T);
end

% Pre-allocate output
X = zeros(1, floor(DUR*SR));

[noscs, ncols] = size(F);
% Implicit zero point at t=0
lastTix = 0;
lastF = F(:,1);
lastM = zeros(noscs,1);
lastP = zeros(noscs, 1);

% Synthesize by columns
for c = 1:ncols
  % Generate all the sample points up to the new time
  Tix = round(T(c)*SR);
  thisM = M(:, c);
  thisF = F(:, c);
  if (Tix > lastTix)
    % Do have some new samples to emit
    npts = Tix - lastTix;
    ix = 0:(npts - 1); % index of each sample in this Tick
    sampts = (lastTix+ix)/SR; % time of each sample from start
    % Make sure fade-outs are to the same frq, not to 0 frq
    thisF = thisF + (thisM == 0 & thisF == 0).*lastF;
	% if current timeslice F and M is zero, add lastFreq to thisfreq
	
	% add to freq & M depending on proportion of distance thru timeslice
    Ft = lastF*ones(1,npts) + (thisF - lastF)*(ix/npts);
    Mt = lastM*ones(1,npts) + (thisM - lastM)*(ix/npts);
	
    Pt = cumsum([lastP, 2*pi*Ft/SR]')';
    lastP = Pt(:, npts+1);
    X(lastTix+1+ix) = sum(Mt.*cos(Pt(:, 1:npts)),1); %einar added ,1)
  end
  lastTix = Tix;
  lastF = thisF;
  lastM = thisM;
end

X=X ./size(F,1);	%Normalize amplitude by number of oscillators
