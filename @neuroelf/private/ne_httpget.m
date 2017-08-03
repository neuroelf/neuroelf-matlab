function varargout = ne_httpget(varargin)
% ne_httpget  - HTTP GET from URL and store to file
%
% FORMAT:       [cont, length = ] ne_httpget(SRC, EVT, URL, ofile, opts)
%
% Input fields:
%
%       SRC, EVT    Matlab-internal handle and keyboard event information
%       URL         URL
%       ofile       output file (or 'var' to retrieve in variable)
%       opts        optional fields
%        .range     1x2 byte range, translated into HTTP Range: bytes=%d-%d
%        .steps     1x2 byte steps, e.g. [1024, length] (Inf is lookup)
%
% Output fields:
%
%       cont        contents of resource (or error)
%       length      length of contents (only valid for range/steps)
%
% Note: this file is built on urlread and urlreadwrite.
%
% Note: if steps is used, several HTTP requests will be made, using the
%       GUI's progress bar to show the download progress, useful for
%       larger files; if the second value is set to Inf, the length
%       will be requested from the HTTP server first, only works well for
%       static documents/files!
%
% Example:
%
%       tmat = ne_httpget(0, 0, ...
%           'http://neuroelf.net/spm8_dartel_template.mat', ...
%           'var', struct('steps', [1048576, Inf]));

%   Matthew J. Simoneau, 13-Nov-2001
%   Copyright 1984-2006 The MathWorks, Inc.
%   $Revision: 1.3.2.8 $ $Date: 2006/12/20 07:16:57 $

% Version:  v1.1
% Build:    16060222
% Date:     Jun-02 2016, 10:33 PM EST
% URL/Info: http://neuroelf.net/

% global access
global ne_gcfg;
try
    ch = ne_gcfg.h;
    shprog = isxfigure(ch.Progress, true);
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    shprog = false;
end

% initialize output
varargout = cell(1, nargout);

% argument check
if nargin < 4 || ~ischar(varargin{3}) || numel(varargin{3}) < 10 || ...
   (~strcmpi(lsqueeze(varargin{3}(1:7))', 'http://') && ...
    ~strcmpi(lsqueeze(varargin{3}(1:8))', 'https://')) || ...
   ~ischar(varargin{4}) || isempty(varargin{4})
    return;
end
url = varargin{3}(:)';
ofile = varargin{4}(:)';
if nargin < 5 || ...
   ~isstruct(varargin{5}) || ...
    numel(varargin{5}) ~= 1
    opts = struct;
else
    opts = varargin{5};
end
if ~isfield(opts, 'range') || ...
   ~isa(opts.range, 'double') || ...
    numel(opts.range) ~= 2 || ...
    any(isinf(opts.range) | isnan(opts.range) | opts.range < 0 | opts.range ~= fix(opts.range)) || ...
    opts.range(1) > opts.range(2)
    opts.range = [];
end
if ~isfield(opts, 'steps') || ...
   ~isa(opts.steps, 'double') || ...
    numel(opts.steps) ~= 2 || ...
    any(isnan(opts.steps) | opts.steps <= 0 | opts.steps ~= fix(opts.steps)) || ...
    isinf(opts.steps(1)) || ...
    opts.steps(1) >= opts.steps(2)
    opts.steps = [];
elseif isinf(opts.steps(2))
    [nullcont, clength] = ne_httpget(0, 0, url, 'var', struct('range', [0, 1]));
    opts.steps(2) = clength;
end

% stepwise
if ~isempty(opts.steps)
    cstep = opts.steps(1);
    clength = opts.steps(2);
    opts.steps = [];
    sts = 0:cstep:clength;
    if sts(end) ~= clength
        sts(end+1) = clength;
    end

    % show progress bar
    cprog = ne_progress(0, 0, ...
        {true, 0, sprintf('Downloading file (0/%dbytes...)', clength)});
    drawnow;

    % loop
    ccont = cell(numel(sts) - 1, 1);
    for cnum = 1:numel(ccont)
        try
            opts.range = sts(cnum:cnum+1);
            opts.range(2) = opts.range(2) - 1;
            ccont{cnum} = ne_httpget(0, 0, url, 'var', opts);
            if numel(ccont{cnum}) ~= (sts(cnum+1) - sts(cnum))
                error( ...
                    'neuroelf:HTTPError', ...
                    'Error retrieving chunk with range (%d-%d) from %s.', ...
                    opts.range(1), opts.range(2), url ...
                );
            end
            if shprog
                ch.Progress.Progress(cnum/numel(ccont), sprintf( ...
                    'Downloading file (%d/%dbytes...)', sts(cnum+1), clength));
            end
        catch ne_eo;
            break;
        end
    end
    
    % deal with bad content
    if any(lsqueeze(cellfun('prodofsize', ccont)) ~= diff(sts(:)))
        if exist('ne_eo', 'var') == 0
            ne_eo = struct('message', 'Incomplete chunks.');
        end
        uiwait(warndlg(sprintf('Error downloading file:\n%s.', ne_eo.message), ...
            'NeuroElf - warning', 'modal'));

    % everything went well
    else
        ccont = lsqueeze(lsqueezecells(ccont));
        ccont = cat(1, ccont{:});
        varargout{1} = ccont;
        varargout{2} = clength;
        if ~strcmp(ofile, 'var')
            binwrite(ofile, ccont);
        end
    end

    % re-set progress bar
    ne_progress(0, 0, cprog);
    return;
end

% this function requires Java.
if ~usejava('jvm')
   error( ...
       'neuroelf:NoJVM', ...
       'Requires Java VM to be loaded.' ...
   );
end

% import required classes
import com.mathworks.mlwidgets.io.InterruptibleStreamCopier;

% be sure the proxy settings are set
com.mathworks.mlwidgets.html.HTMLPrefs.setProxySettings

% create a urlConnection
[urlConnection, errorid, errormsg] = urlreadwrite(url, opts);
if isempty(urlConnection)
    error( ...
        'neuroelf:JavaError', ...
        'Error creating connection: (%d, %s).', ...
        errorid, errormsg ...
    );
end

% second output, try to access all length
if nargout > 1
    crange = urlConnection.getHeaderField('Content-Range');
    if ~isempty(crange)
        varargout{2} = str2double(regexprep(char(crange), '^.*\/', ''));
    else
        varargout{2} = -1;
    end
end

% read the data from the connection
try
    inputStream = urlConnection.getInputStream;

    % copy to output stream
    byteArrayOutputStream = java.io.ByteArrayOutputStream;
    isc = InterruptibleStreamCopier.getInterruptibleStreamCopier;
    isc.copyStream(inputStream,byteArrayOutputStream);

    % close streams
    inputStream.close;
    byteArrayOutputStream.close;

    % copy to output
    output = typecast(byteArrayOutputStream.toByteArray', 'uint8');
catch ne_eo;
    rethrow(ne_eo);
end

% make output public
if nargout > 0
    varargout{1} = output(:);
end

% write to file
try
    if ~strcmp(ofile, 'var')
        binwrite(ofile, output);
    end
catch ne_eo;
    rethrow(ne_eo);
end


function [urlConnection, errorid, errormsg] = urlreadwrite(urlChar, opts)
urlConnection = [];
errorid = '';
errormsg = '';
try
    handler = sun.net.www.protocol.http.Handler;
catch exception %#ok
    handler = [];
end
try
    if isempty(handler)
        url = java.net.URL(urlChar);
    else
        url = java.net.URL([],urlChar,handler);
    end
catch exception %#ok
    errorid = 'MATLAB:urlread:InvalidUrl';
    errormsg = 'Either this URL could not be parsed or the protocol is not supported.';
    return
end
mwtcp = com.mathworks.net.transport.MWTransportClientPropertiesFactory.create();
proxy = mwtcp.getProxy();
if isempty(proxy)
    urlConnection = url.openConnection;
else
    urlConnection = url.openConnection(proxy);
end
if ~isempty(opts.range)
    urlConnection.setDoOutput(true);
    urlConnection.setRequestProperty( ...
        'Range',sprintf('bytes=%d-%d', opts.range(1), opts.range(2)));
end
