% FUNCTION ne_setwindowpos: get new position from text control
function varargout = ne_setwindowpos(varargin)

% Version:  v1.1
% Build:    16061321
% Date:     Jun-13 2016, 9:40 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2015, 2016, Jochen Weber
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

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% get list of text objects (Matlab handles)
if nargin < 3 || ...
   ~ischar(varargin{3}) || ...
    isempty(varargin{3}) || ...
   ~isfield(ne_gcfg.cc, varargin{3}(:)')
    fromroot = true;
    cc = ne_gcfg.fcfg;
    ch = ne_gcfg.h;
else
    fromroot = false;
    iSat = varargin{3}(:)';
    ch = ne_gcfg.cc.(iSat);
    cc = ch.Config;
end
tobj = cc.txtpos;

% if the callback object is not one of those
if nargin > 0 && ...
    numel(varargin{1}) == 1 && ...
    ishandle(varargin{1})
    rcbo = varargin{1};
else
    rcbo = gcbo;
end
cobj = find(tobj == rcbo);
if isempty(cobj)

    % return already
    return;
end

% get current (former) position
cpos = cc.cpos;

% get new position (one of the three ordinates)
cstr = get(tobj(cobj), 'String');

% make check depending on control, and if bad set old value -> TAL
if cobj <= 3 && ...
   (isempty(regexpi(cstr, '^[\+\-]?\d+(\.\d)?$')) ||  ...
    abs(str2double(cstr(:)')) > 128)
    cstr = sprintf('%.0f', cpos(cobj));
    set(tobj(cobj), 'String', cstr);

% -> BVS
elseif cobj > 3 && ...
   (isempty(regexpi(cstr, '^\+?\d+(\.\d)?$')) ||  ...
    str2double(cstr(:)') <= 0 || ...
    str2double(cstr(:)') > 512)

    % we need the slicing var
    svar = cc.SliceVar;

    % for no object, assume VMR logic
    if ~isxff(svar, true)

        % 129 - value
        cstr = sprintf('%.0f', 129 - cpos(cobj - 3));

    % for actual object
    else

        % convert to voxel coordinates
        bvpos = svar.RunTimeVars.Trf' * [cpos(:); 1];
        if bvpos(cobj - 3) == fix(bvpos(cobj - 3))
            cstr = sprintf('%d', bvpos(cobj - 3));
        else
            cstr = sprintf('%g', 0.1 * round(10 * bvpos(cobj - 3)));
        end
    end
    set(tobj(cobj), 'String', cstr);

% for good values
else

    % depending on object
    switch (cobj)

        % TAL coordinates
        case {1, 2, 3}

            % as is
            cpos(cobj) = str2double(cstr(:)');

        % BVS coordinates
        case {4, 5, 6}

            % we need the slicing object
            svar = cc.SliceVar;

            % for no object, assume VMR logic
            if ~isxff(svar, true)

                % 128 - value
                cpos(cobj - 3) = 129 - str2double(cstr(:)');
            else

                % get the number of all three objects
                try
                    dpos = [str2double(get(ch.Coord.VEdX, 'String')), ...
                            str2double(get(ch.Coord.VEdY, 'String')), ...
                            str2double(get(ch.Coord.VEdZ, 'String'))];
                    if numel(dpos) ~= 3 || ...
                        any(isinf(dpos) | isnan(dpos) | dpos < -256 | dpos > max(svar.GetVolumeSize))
                        error('BAD_DPOS');
                    end

                    % now do the back transformation
                    if ~isequal(cc.strans, eye(4))
                        trans = cc.strans;
                        if cc.orient == 'n'
                           trans(:, 1) = -trans(:,1);
                        end
                    else
                        if cc.orient ~= 'n'
                            trans = [];
                        else
                            trans = eye(4);
                            trans(1) = -1;
                        end
                    end

                    % use Trf to get voxel position
                    tplus = svar.RunTimeVars.TrfPlus';
                    itrf = inv(svar.RunTimeVars.Trf')';
                    if ~isempty(trans)
                        cpos = [dpos, 1] * itrf * trans' * tplus;
                    else
                        cpos = [dpos, 1] * itrf * tplus;
                    end
                    cpos(end) = [];

                    % type check
                    if isxff(svar, 'vmr')
                        if any([svar.VoxResX, svar.VoxResY, svar.VoxResZ] > 0.5)
                            cpos = round(cpos);
                        else
                            cpos = 0.5 .* round(2 .* cpos);
                        end
                    else
                        cpos = 0.1 .* round(10 .* cpos);
                    end

                catch ne_eo;
                    ne_gcfg.c.lasterr = ne_eo;
                end
                if cc.orient == 'n'
                    cpos(1) = -cpos(1);
                end
            end
        otherwise
            return;
    end

    % linked browsing
    if ne_gcfg.c.linked

        % set in all windows
        ne_setslicepos(0, 0, cpos, 'OnPosition');

    % from root (main ) UI
    elseif fromroot

        % set back in figure configuration
        ne_gcfg.fcfg.cpos = cpos;

        % update window
        ne_setslicepos(0, 0, [], 'OnPosition');

    % from satellite
    else

        % set back in satellite configuration
        ne_gcfg.cc.(iSat).Config.cpos = cpos;

        % update window
        ne_setsatslicepos(0, 0, iSat);
    end
end
