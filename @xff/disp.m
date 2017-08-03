function disp(xo)
%DISP  Display method for XFF objects.
%   DISP(OBJ) displays the content of an XFF object. For the ROOT object
%   (factory), a list of files is being displayed.
%
%   See also XFF

% Version:  v1.1
% Build:    16012414
% Date:     Jan-24 2016, 2:50 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2014, 2016, Jochen Weber
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

% global storage and require methods
global xffsngl ne_methods;

% for non-singleton objects
if numel(xo) ~= 1

    % just display typical size information
    osz = sprintf('%d-by-', size(xo));
    fprintf('  %s <a href="matlab:help xff">xff</a> array.\n\n', osz(1:end-4));

% invalid object
elseif ~isvalid(xo)
    fprintf('  handle to deleted xff\n\n');

% valid object
else
    try

        % object banner
        if xo.L(1) == 'X'
            disp('  <a href="matlab:help xff">xff</a> factory contents:');
            disp(' ');
        else
            otype = upper(xo.S.Extensions{1});
            disp(['  <a href="matlab:help xff">xff</a> (' otype ') object with properties:']);
            disp(' ');
        end

        % content
        ne_methods.dispstruct(xo.C);

        % if the current object is the ROOT/factory
        if xo.L(1) == 'X'

            % how much memory does the singleton require
            tspc = whos('xffsngl');
            tspc = tspc.bytes;

            % display some information about loaded files
            disp(' ');
            disp('  List of currently loaded objects:');
            disp(repmat('-', 1, 80));
            disp('   # | Type  | Mem (MB) | Filename or [Document ID]');
            disp(repmat('-', 1, 80));
            for oc = 2:size(xffsngl.OBJS, 1)
                tspo = xffsngl.OBJS{oc, 4};
                tspoc = tspo.C;
                tsps = whos('tspoc');
                tspoc = sprintf('%8.3f', tsps.bytes / 1048576);

                % without filename
                if isempty(tspo.F)
                    fprintf('%4d | %4s  | %s | [%s]\n', oc - 1, tspo.S.Extensions{1}, tspoc(1:8), tspo.L);

                % with filename
                else
                    fprintf('%4d | %4s  | %s | %s\n', oc - 1, tspo.S.Extensions{1}, tspoc(1:8), tspo.F);
                end

                % memory allocation
                tspc = tspc + tsps.bytes;
            end

            % print total memory consumption
            if size(xffsngl.OBJS, 1) > 1
                disp(repmat('-', 1, 80));
            end
            fprintf('  Total memory occupied: %.3f MByte\n', tspc / 1048576);
        end
        disp(' ');

    % deal with errors
    catch xfferror
        rethrow(xfferror);
    end
end
