function [cpstatus,cpmsg] = cpfile(fromfile,tofile,force)
% cpfile  - copies a file
%
% FORMAT:       [cpstatus,cpmsg] = cpfile(fromfile, tofile [,force])
%
% Input Fields:
%
%       fromfile    name of the file to use as source file
%       tofile      name of the target file
%       force       if given and set to 1, 'force', or 'overwrite'
%                   copying will be forced
%
% if requested, returns a PERL-like status (0 for unsuccessful!) and
% a message indicating what went wrong (or empty if successful)

% Version:  v0.9b
% Build:    11050712
% Date:     Apr-09 2011, 1:55 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, Jochen Weber
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%     * Neither the name of Columbia University nor the
%       names of its contributors may be used to endorse or promote products
%       derived from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% settings
icpbuffer = 65536;
cpstatus  = 0;
cpmsg     = '';

% enough arguments ?
if nargin < 2
    error( ...
        'neuroelf:TooFewArguments',...
        'Too few arguments. Try ''help %s''.',...
        mfilename ...
    );
end

if ~ischar(fromfile) || ...
   ~ischar(tofile) || ...
    isempty(fromfile) || ...
    isempty(tofile)
    error( ...
        'neuroelf:BadArgument',...
        'Bad argument.' ...
    );
end

% options
if nargin < 3 || ...
    isempty(force)
    force = 0;
elseif isnumeric(force) && ...
   ~isnan(force(1)) && ...
    force(1) ~= 0
    force = 1;
elseif ischar(force) && ...
    any('fo' == lower(force(1)))
    force = 1;
else
    force = 0;
end

% filename mangling check
if ispc
    fromfile = strrep(fromfile(:)', '/', filesep);
    tofile = strrep(tofile(:)', '/', filesep);
else
    fromfile = strrep(fromfile(:)', '\', filesep);
    tofile = strrep(tofile(:)', '\', filesep);
end

try

    % check for target being a folder
    if exist(tofile,'dir') == 7
        [topartfn{1:3}] = fileparts(fromfile);
        if ~strcmp(filesep, tofile((end+1-length(filesep)):end))
            tofile = [tofile filesep];
        end
        tofile = [tofile topartfn{2} topartfn{3}];
    end

    % check for existing target file
    if exist(tofile, 'file') == 2 && ...
        force < 1
        error( ...
            'neuroelf:NoForceRequested',...
            'Tofile ''%s'' already exists; not overwritten.',...
            tofile ...
        );
    end

    % open source file and check for handle
    sfp = fopen(fromfile, 'r');
    if sfp < 1
        error( ...
            'neuroelf:FileNotReadable',...
            'Couldn''t read from file: %s',...
            fromfile ...
        );
    end

    % open destination file and check for handle
    tfp = fopen(tofile, 'w');
    if tfp < 1
        fclose(sfp);
        error( ...
            'neuroelf:FileNotWritable',...
            'Couldn''t write to file: %s',...
            tofile ...
        );
    end

    % rewind source file, write content to target and close file handles
    frewind(sfp);
    while ~feof(sfp)
        buff = fread(sfp, icpbuffer, '*uint8');
        fwrite(tfp, buff, 'uint8');
        if numel(buff) < icpbuffer
            break;
        end
    end
    fclose(sfp);
    fclose(tfp);

    cpstatus = 1;

catch ne_eo;
    cpmsg = ne_eo.message;
end
