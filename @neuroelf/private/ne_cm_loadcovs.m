% FUNCTION ne_cm_loadcovs: load covariates
function ne_cm_loadcovs(varargin)

% Version:  v0.9b
% Build:    11051315
% Date:     Apr-09 2011, 11:08 PM EST
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

% global variable
global ne_gcfg;
cc = ne_gcfg.fcfg.CM;
ch = ne_gcfg.h.CM.h;

% request file
[cvfile, cvpath] = uigetfile( ...
    {'*.cov',       'COVariate files (*.cov)'; ...
     '*.rtc;*.sdm', 'Design matrix files (*.rtc, *.sdm)'; ...
     '*.mat',       'MAT-files (*.mat)'; ...
     '*.txt',       'Text files (*.txt)'; ...
     '*.*',         'All files (*.*)'}, ...
     'Please select the file containing the covariates...', ...
     'MultiSelect', 'off');

% file selected
if isequal(cvfile, 0) || ...
    isequal(cvpath, 0) || ...
    numel(cvfile) < 5
    return;
end
cvfull = [cvpath, cvfile];

% try to load file
ext = lower(cvfile(end-2:end));
try

    % depending on extension
    switch (ext)

        % for COV files
        case {'cov'}

            % use xff
            cvobj = xff(cvfull);

            % and get covariates
            cvcont = cvobj.Covariates;
            cvobj.ClearObject;

            % then create names
            cvname = cell(1, size(cvcont, 2));
            for cvc = 1:numel(cvname)
                cvname{cvc} = sprintf('%s - Covariate %d', cvfile, cvc);
            end

        % for RTC/SDM files
        case {'rtc', 'sdm'}

            % follow the same logic
            cvobj = xff(cvfull);
            cvcont = cvobj.RTCMatrix;
            cvname = cvobj.PredictorNames(1:size(cvcont, 2));
            cvobj.ClearObject;

        % for all other files
        otherwise

            % simply try to use "load" (will fail and be caught...)
            cvobj = load(cvfull);

            % for numeric variables
            if isnumeric(cvobj)

                % simply take numbers
                cvcont = double(cvobj);

                % and create names
                cvname = cell(1, size(cvcont, 2));
                for cvc = 1:numel(cvname)
                    cvname{cvc} = sprintf('%s - Covariate %d', cvfile, cvc);
                end

            % for a structure (MAT-file with only one variable)
            elseif isstruct(cvobj) && ...
                numel(fieldnames(cvobj)) == 1

                % get the variable name
                cvofld = fieldnames(cvobj);

                % and content
                cvcont = cvobj.(cvofld{1});

                % which must be numeric
                if ~isnumeric(cvcont)
                    error('Unsupported MAT file content.');
                end
                cvcont = double(cvcont);

                % and then create names
                cvname = cell(1, size(cvcont, 2));
                for cvc = 1:numel(cvname)
                    cvname{cvc} = sprintf('%s - %s %d', cvfile, cvofld{1}, cvc);
                end

            % for a structure (MAT-file with exactly two variables)
            elseif isstruct(cvobj) && ...
                numel(fieldnames(cvobj)) == 2

                % get the variable name
                cvofld = fieldnames(cvobj);

                % and content
                cvcont1 = cvobj.(cvofld{1});
                cvcont2 = cvobj.(cvofld{2});
                if ~isnumeric(cvcont1)
                    cvcont2 = cvcont1;
                    cvcont1 = cvobj.(cvofld{2});
                end

                % which must be numeric
                if ~isnumeric(cvcont1) || ...
                   ~iscell(cvcont2) || ...
                    numel(cvcont2) ~= size(cvcont1, 2)
                    error('Unsupported MAT file content.');
                end
                cvcont = double(cvcont1);

                % and then check/create names
                cvname = cvcont2(:)';
                for cvc = 1:numel(cvname)
                    if ~ischar(cvname{cvc}) || ...
                        isempty(cvname{cvc}) || ...
                       (cvc > 1 && ...
                        any(strcmpi(cvname(1:cvc-1), cvname{cvc})))
                        cvname{cvc} = sprintf('%s - %s %d', cvfile, cvofld{1}, cvc);
                    end
                end

            % other forms of content not supported
            else
                error('Unsupported MAT/TXT file encountered.');
            end
    end

% inform the user in case of problems
catch ne_eo;
    uiwait(msgbox(sprintf('Error loading the covariate file: %s', ne_eo.message), ...
        'NeuroElf GUI - info', 'modal'));
    return;
end

% incorrect size?
if size(cvcont, 1) ~= cc.nsubs

    % let user name
    uiwait(msgbox('The covariate has the wrong number of values.', ...
        'NeuroElf GUI - info', 'modal'));

    % and don't do anything
    return;
end

% otherwise add to list
cvc = size(cc.covs, 1);
cc.covs(cvc+1:cvc+numel(cvname), 1) = cvname(:);
for ncvc = 1:numel(cvname)
    cc.covs{cvc + ncvc, 2} = cvcont(:, ncvc);
end

% then set in handle
ch.Covs.String = cc.covs(:, 1);

% and update global array
ne_gcfg.fcfg.CM = cc;

% update in GLM
ne_cm_updatertv;
