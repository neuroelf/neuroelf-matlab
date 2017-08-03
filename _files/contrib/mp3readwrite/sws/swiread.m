function [F,M,timeslice, nOscs] = swiread(FILENAME)

%	 Read a Haskins-format sinewave speechdata file.
%	FILENAME is the name of a text data file containing the frequencyand magnitude
%	parameters for sinewave synthesis.  Result is F and M matrices suitable for synthtrax.m.
%	Timeslice is the time index  for regular and irregularly sampled files, suitable for synthtrax.
%
%	dpwe@icsi.berkeley.edu 1996aug22 Revised: 1998Sept18
%	Also revised by Steve Frost 9.98 and 1.00
%
% SWI files have the format:
%
%    Number of oscillators
%      Time0
%         frq,mag   for 1st oscillator
%         frq,mag   for 2nd oscillator
%         .. for as many oscillators as specified
%      Time1
%         frq,mag  ... etc.
% Times are in ms, frq in Hz, mag in linear units
%
%% SWX files have the format:
%
%    Number of oscillators
%      Time0	oscillator1freq	oscillator1mag	oscillator2freq	oscillator2mag.....oscillatorNfreq	oscillatorNmag
%      Time1	oscillator1freq	oscillator1mag	oscillator2freq	oscillator2mag.....oscillatorNfreq	oscillatorNmag
%         for as many Time values specified
% Times are in ms, frq in Hz, mag in linear units
%



%	Process data for synthtrax
colchunk = 100;
col = 0;

fid = fopen(FILENAME, 'r');
if (fid == -1)
  fprintf(1, 'readswi: unable to read %s\n', FILENAME);
else
  nOscs = fscanf(fid, '%d', 1);
% 	Increase the arrays in chunks of colchunk cols to avoid slow
% 	matrix growing.
  emptyF = zeros(nOscs, colchunk);
  F = emptyF;
  M = emptyF;
  Fcols = colchunk;
  T = emptyF(1,:);

  endoffile = 0;
  while (endoffile == 0)
    [time,count] = fscanf(fid, '%f', 1);

    if (count == 0)
      endoffile = 1;
    else
      col = col+1;
      if(col > Fcols)
        % We ran out of empty columns - grow the matrices
        F = [F, emptyF];
        M = [M, emptyF];
        T = [T, emptyF(1,:)];
        Fcols = Fcols + colchunk;
      end
      for osc = 1:nOscs
        F(osc,col) = fscanf(fid, '%f,', 1);
        M(osc,col) = fscanf(fid, '%f',1);
      end
      T(col) = time;
    end
  end
  fclose(fid);
  
  % Trim off excess empty columns
  F = F(:,1:col);
  M = M(:,1:col);
  T = T(1:col);
  timeslice=T';
end

