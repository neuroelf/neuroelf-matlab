function [varargout] = neuroelf(varargin)
% neuroelf (Object Class)
%
% The neuroelf object is a singleton object (factory-style, a bit like the
% JavaScript Math object/class) that provides a hub for the functionality
% in the NeuroElf toolbox, such that
%
% n = neuroelf;
% f = n.findfiles(ARGUMENTS{:});
%
% does the same as in previous versions
%
% f = findfiles(ARGUMENTS{:});
%
% would have done. Additional methods and properties will be added in time.
%
% Overhead for this kind of "hidden implementation" is approximately 1e-4s
% per call to a method; so this syntax is to be considered for relatively
% high-level methods (functions) that are not called repeatedly in a loop.
% As an alternative, you can request a function handle using
%
% findfiles = n.findfiles;
% f = findfiles(ARGUMENTS{:});
%
% which then has an overhead approximately < 1-e6s per call.
%
% To establish several function handles (and without assigning an
% additional variable), the following syntax assigns method names to the
% current (caller) workspace:
%
% using(neuroelf, {'findfiles', 'mfileparts', 'multimatch'});
%
% such that the variables findfiles, mfileparts, and multimatch hold the
% function handles to those three functions.
%
% Finally, all methods are also available via the global ne_methods struct:
%
% n = neuroelf;
% global ne_methods;
% findfiles = ne_methods.findfiles;
% f = findfiles(ARGUMENTS{:});
%
% To access the help of any of the functions, use
%
% n = neuroelf;
% n.help('findfiles');
%
% To allow access to function also without creating an additional (object)
% variable, function names (as well as arguments) can also be passed
% directly into the class constructor function (in which case the return
% value will not be of the class type!). Example:
%
% f = neuroelf('findfiles', pwd, '*.ext');

% Version:  v1.1
% Build:    16053016
% Date:     May-30 2016, 4:29 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2014 - 2016, Jochen Weber
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

% global methods
global ne_methods;

% persistent internal storage for factory
persistent neos;

% initialize (once)
if numel(neos) ~= 1 || ~isstruct(neos)

    % set to struct
    neos = struct('cfg', struct, 'init', false, 'meth', struct, 'mets', struct);

    % find methods (manually adding findfiles)
    npath = neuroelf_path;
    mpath = [npath filesep '@neuroelf' filesep 'private'];
    meths = sort(cat(1, ...
        lsqueeze(findfiles(mpath, '*.m', 'depth=1', 'relative=')), ...
        lsqueeze(findfiles(mpath, ['.' lower(mexext)], 'depth=1', 'relative='))));

    % for each method
    for mc = 1:numel(meths)

        % get shortname (without extension)
        smeth = regexprep(meths{mc}, ['\.(m|' lower(mexext) ')$'], '');

        % if an "M-file"
        if strcmp(meths{mc}(end-1:end), '.m')

            % read content
            mcont = asciiread([mpath filesep meths{mc}]);

            % replace old line delimiter (Windows) to char(10) only
            if ~isempty(strfind(mcont, char([13, 10])))
                mcont(mcont == 13) = [];
            end

            % split to lines (for arguments and help)
            mcont = splittocellc(mcont, char(10));

            % remove "^\s*function\s+" from header (first line)
            mhead = regexprep(mcont{1}, '^\s*function\s+', '');

            % any assignment (outputs)?
            if any(mhead == '=')

                % only look at that part
                mouts = ddeblank(regexprep(mhead, '^([^=]+)\s*=.*$', '$1'));

                % potentially multiple outputs?
                if mouts(1) == '[' && mouts(end) == ']'

                    % don't worry about [ and ] though!
                    mouts = mouts(2:end-1);
                end

                % simply split by , and ' ' (as many as possible)
                mouts = splittocellc(mouts, ', ', true, true);

            % no outputs
            else
                mouts = {};
            end

            % number of outputs
            nmouts = numel(mouts);

            % extended to Inf, if varargout!
            if nmouts > 0 && strcmp(mouts{end}, 'varargout')
                nmouts = Inf;
            end

            % inputs?
            if any(mhead == '(')

                % only look at that part
                mins = ddeblank(regexprep(mhead, '^[^\(]+\(([^\)]+)\)\s*$', '$1'));

                % but make sure it's valid!
                if any(mins == '(' | mins == ')')
                    mins = {};
                else

                    % split to components (like outputs)
                    mins = splittocellc(mins, ', ', true, true);
                end
            else
                mins = {};
            end

            % number of inputs
            nmins = numel(mins);

            % extended to Inf, if varargin!
            if nmins > 0 && strcmp(mins{end}, 'varargin')
                nmins = Inf;
            end

            % get lines with comments (help candidates)
            mhelp = find(~cellfun('isempty', regexp(mcont, '^\%')));

            % some help found!
            if ~isempty(mhelp)

                % remove lines after first block
                mhelp((1+findfirst(diff(mhelp) > 1)):end) = [];
            end

            % if any help (left)
            if ~isempty(mhelp)

                % store as single string (glue back together, remove '% ')
                mhelp = gluetostringc(regexprep(mcont(mhelp), '^\%\s?', ''), char(10), true);

            % without help
            else

                % construct list of inputs
                if nmins > 0
                    sins = sprintf('%s, ', mins{:});
                    sins = ['(' sins(1:end-2) ')'];
                else
                    sins = '';
                end

                % construct list of outputs
                if nmouts > 0
                    if nmouts == 1
                        souts = [mouts{1} ' = '];
                    else
                        souts = sprintf('%s, ', mouts{:});
                        souts = ['[' souts(1:end-2) '] = '];
                    end
                else
                    souts = '';
                end

                % construct call (prototypically)
                mhelp = sprintf('%s%s%s', souts, smeth, sins);
            end

            % store in methods array
            neos.meth.(smeth) = {eval(['@' smeth]), nmouts, mouts, nmins, mins, mhelp};

        % for MEX files
        else

            % if no M-file with same name found
            if ~isfield(neos.meth, smeth)

                % assume infinite inputs and outputs
                neos.meth.(smeth) = {eval(['@' smeth]), Inf, {'varargout'}, ...
                    Inf, {'varargin'}, sprintf('[varargout] = %s(varargin)', smeth)};
            end
        end

        % store shortcut
        neos.mets.(smeth) = neos.meth.(smeth){1};
    end

    % initialize
    neos.init = true;

    % then store methods in global variable as well
    ne_methods = neos.mets;

    % finally, if the toolbox hasn't been setup yet, run setup
    if any(npath(end) == '0123456789')
        npbuild = regexprep(npath, '^.*_(\d+)$', '$1');
    else
        npbuild = 'release';
    end
    sdone = [neuroelf_path('cache') filesep 'setup_' npbuild '.done'];
    if numel(dbstack) == 1 && ...
        (exist(sdone, 'file') ~= 2 || ~strcmp(ddeblank(asciiread(sdone)), 'done'))
        neuroelf_setup;
    end
end

% if outputs requested or a function call is intended
if nargout > 0 || (nargin > 0 && ischar(varargin{1}) && isfield(neos.mets, varargin{1}(:)'))

    % preset outputs
    varargout = cell(1, max(1, nargout));

    % return the singleton object if no inputs given (or not a function)
    if nargin == 0 || ~ischar(varargin{1}) || isempty(varargin{1}) || ...
       ~isfield(neos.mets, varargin{1}(:)')
        varargout{1} = class(neos, 'neuroelf');
        if nargout > 1
            varargout{2} = neos.mets;
        end

    % otherwise
    else

        % get the method's configuration
        meth = neos.meth.(varargin{1}(:)');

        % allow for errors
        try

            % if the method has outputs and outputs are requested
            if meth{2} > 0 && nargout > 0

                % call with outputs
                [varargout{1:min(meth{2}, nargout)}] = ...
                    feval(neos.mets.(varargin{1}(:)'), varargin{2:end});

            % otherwise
            else

                % call without outputs
                feval(neos.mets.(varargin{1}(:)'), varargin{2:end});
            end

        % deal with errors
        catch ne_eo;
            rethrow(ne_eo);
        end
    end

% without and inputs AND outputs
elseif nargin == 0

    % open GUI
    neuroelf_gui;
end
