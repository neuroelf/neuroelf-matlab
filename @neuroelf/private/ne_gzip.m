function [varargout] = ne_gzip(varargin)
% ne_gzip  - console/shell interface for MATLAB
%
% a simple interface to gzip functionality that *should* work under
% MS-Windows as well (please download and install the binary
% e.g. courtesy of www.gzip.org).
%
% since all arguments must be strings, you can use both command and
% functional calling convention in MATLAB. as a function, you can get
% both return status and console output, only return status, or nothing.
%
% ne_gzip('arg1','arg2',...);
% status = ne_gzip('arguments');
% [status,output] = ne_gzip('arg1','arg2',...);
%
% all arguments are glued together with a blank, so it doesn't matter
% whether you call with one argument or several.
%
% examples:
%
% ne_gzip -9 /tmp/hallo.bmp
% s = ne_gzip('-d /home/user123/myfiles/vol1.img.gz');
% ne_gzip -t c:\temp\file123.gz
% for i = 1:nfiles, ne_gzip('-9',files{i}); end
%
% when nothing is passed out of this interface, the text output will be
% shown (disp-ed), but no return value will be available; this also
% applies to the command-style call. with one or more return values
% nothing will be printed to the MATLAB command window.
%
% no input arguments do a call to gzip's help page (gzip -h).

% Version:  v0.9c
% Build:    11050712
% Date:     May-02 2011, 5:03 PM EST
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

% default argument
callargs=' -h';

% argument given
if nargin > 0

    % argument is not char
    if ~ischar(varargin{1})

        % only valid if true output is requested and UINT8 input
        if nargout < 1 || ...
            nargin < 2 || ...
           ~isa(varargin{1}, 'uint8')
            error( ...
                'neuroelf:BadArgument', ...
                'Only uint8 content can be handled via this syntax.' ...
            );
        end

        % also check option argument
        if ~ischar(varargin{2}) || ...
            numel(varargin{2}) ~= 2 || ...
            varargin{2}(1) ~= '-' || ...
           ~any('d123456789' == lower(varargin{2}(2)))
            error( ...
                'neuroelf:BadArgument', ...
                'Illegal gzip mode for stream (de)compression.' ...
            );
        end

        % write content to file
        tn = tempname;

        % with correct extension (-1...9 / -d)
        if numel(varargin{2}) == 2 && ...
            all(lower(varargin{2}) == '-d')
            itn = [tn '.gz'];
            otn = tn;
        else
            itn = tn;
            otn = [tn '.gz'];
        end

        % try opening/writing file
        gfp = fopen(itn, 'w');
        if gfp < 1
            error( ...
                'neuroelf:FileNotWritable', ...
                'Couldn''t write temporary file (%s).', ...
                itn ...
            );
        end
        fwrite(gfp, varargin{1}, 'uint8');
        fclose(gfp);

        % recall gzip with filename and option
        [s, w] = ne_gzip(varargin{2}, itn);

        % check status
        if s ~= 0

            % delete file on failure
            try
                delete(itn);
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
            end

            % bail out now
            error( ...
                'neuroelf:ErrorCallingExecutable', ...
                'Error processing stream (de)compression: %s.', ...
                w ...
            );
        end

        % no failure -> open result file
        gfp = fopen(otn, 'r');

        % on error
        if gfp < 1

            % delete both files and bail out
            try
                delete(otn);
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
            end
            try
                delete(itn);
            catch ne_eo;
                neuroelf_lasterr(ne_eo);
            end
            error( ...
                'neuroelf:FileNotReadable', ...
                'Couldn''t read gzip output file (%s).', ...
                otn ...
            );
        end

        % read, close, and delete file
        varargout{1} = fread(gfp, Inf, '*uint8');
        fclose(gfp);
        try
            delete(otn);
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
        end
        return;
    end

    % first argument must then be char
    callargs = ' ';

    % parse arguments into chair
    for c = 1:nargin

        % only allow char arguments though
        if ~ischar(varargin{c})
            error( ...
                'neuroelf:BadArgument', ...
                'Non-character argument found.' ...
            );
        end

        % put into chain
        callargs = [callargs varargin{c} ' '];
    end
end

% call gzip utility
[s, r] = system(['gzip ' callargs]);

% no output requested
if nargout == 0

    % show any error
    if ~isempty(r)
        disp(r);
    end

% otherwise
else

    % first argument is result
    varargout{1} = s;

    % second is error message
    if nargout > 1
        varargout{2} = r;
    end
end
