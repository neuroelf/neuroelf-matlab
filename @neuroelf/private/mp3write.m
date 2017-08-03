function mp3write(D,SR,NBITS,FILE,OPTIONS)
% MP3WRITE   Write MP3 file by use of external binary
%   MP3WRITE(Y,FS,NBITS,FILE) writes waveform data Y to mp3-encoded
%     file FILE at sampling rate FS using bitdepth NBITS.  
%     The syntax exactly mirrors AUDIOWRITE.  NBITS must be 16.
%   MP3WRITE(Y,FS,FILE) assumes NBITS is 16
%   MP3WRITE(Y,FILE) further assumes FS = 8000.
%
%   MP3WRITE(..., OPTIONS) specifies additional compression control 
%     options as a string passed directly to the lame encoder
%     program; default is '--quiet -h' for high-quality model.
%
%   Example: 
%   To convert a wav file to mp3 (assuming the sample rate is 
%   supported):
%     [Y,FS] = audioread('piano.wav');
%     mp3write(Y,FS,'piano.mp3');
%   To force lame to use 160 kbps (instead of default 128 kbps)
%   with the default filename extension (mp3):
%     mp3write(Y,FS,'piano','--quiet -h -b 160');
%
%   Note: The actual mp3 encoding is done by an external binary, 
%     lame, which is available for multiple platforms.  Usable
%     binaries are available from:
%     http://labrosa.ee.columbia.edu/matlab/mp3read.html
%
%   Note: MP3WRITE will use the mex file popenw, if available, to
%     open a pipe to the lame encoder.  Otherwise, it will have to
%     write a large temporary file, then execute lame on that file.
%     popenw is available at: 
%       http://labrosa.ee.columbia.edu/matlab/popenrw.html
%     This is a nice way to save large audio files as the
%     incremental output of your code, but you'll have to adapt the
%     central loop of this function (rather than using it directly).
%
%   See also: mp3read, audiowrite, popenw.

% 2005-11-10 Original version
% 2007-02-04 Modified to exactly match audiowrite syntax, and to
%   automatically find architecture-dependent binaries.
% 2007-07-26 Writing of stereo files via tmp file fixed (thx Yu-ching Lin)
%
% $Header: /Users/dpwe/matlab/columbiafns/RCS/mp3write.m,v 1.2 2007/07/26 15:09:16 dpwe Exp $

% find our baseline directory
p = [neuroelf_path('contrib') filesep 'mp3readwrite'];

% %%%%% Directory for temporary file (if needed)
% % Try to read from environment, or use /tmp if it exists, or use CWD
tmpdir = getenv('TMPDIR');
if isempty(tmpdir) || exist(tmpdir,'file')==0
  tmpdir = '/tmp';
end
if exist(tmpdir,'file')==0
  tmpdir = '';
end
% ensure it exists
%if length(tmpdir) > 0 && exist(tmpdir,'file')==0
%  mkdir(tmpdir);
%end

%%%%%% Command to delete temporary file (if needed)
rmcmd = 'rm';

%%%%%% Location of the binary - attempt to choose automatically
%%%%%% (or edit to be hard-coded for your installation)
ext = lower(computer);
if ispc
  ext = 'exe';
  rmcmd = 'del';
end
lame = fullfile(p,['lame.',ext]);

%%%% Process input arguments
% Do we have NBITS?
mynargin = nargin;
if ischar(NBITS)
  % NBITS is a string i.e. it's actually the filename
  if mynargin > 3
    OPTIONS = FILE;
  end
  FILE = NBITS;
  NBITS = 16;
  % it's as if NBITS had been specified...
  mynargin = mynargin + 1;
end

if mynargin < 5
  OPTIONS = '--quiet -h';  % -h means high-quality psych model
end

[nr, nc] = size(D);
if nc < nr
  D = D';
  [nr, nc] = size(D);
end
% Now rows are channels, cols are time frames (so interleaving is right)

%%%%% add extension if none (like audioread)
[path,file,ext] = fileparts(FILE);
if isempty(ext)
  FILE = [FILE, '.mp3'];
end

nchan = nr;
nfrm = nc;

if nchan == 1
  monostring = ' -m m';
else
  monostring = '';
end

lameopts = [' ', OPTIONS, monostring, ' '];

tmpfile = fullfile(tmpdir,['tmp',num2str(round(1000*rand(1))),'.wav']);

audiowrite(tmpfile,D',SR);

cmd = ['"',lame,'"', lameopts, '"',tmpfile, '" "', FILE, '"'];

mysystem(cmd);

% Delete tmp file
mysystem([rmcmd, ' "', tmpfile,'"']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w = mysystem(cmd)
% Run system command; report error; strip all but last line
[s,w] = system(cmd);
if s ~= 0 
  error(['unable to execute ',cmd,' (',w,')']);
end
% Keep just final line
w = w((1+max([0,findstr(w,10)])):end);
% Debug
%disp([cmd,' -> ','*',w,'*']);
