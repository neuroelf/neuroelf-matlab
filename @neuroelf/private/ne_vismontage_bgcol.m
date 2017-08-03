% FUNCTION ne_vismontage_bgcol: handle background color stuff
function ne_vismontage_bgcol(varargin)

% Version:  v0.9b
% Build:    11050712
% Date:     Apr-10 2011, 4:52 PM EST
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

% get tags
tags = ne_gcfg.h.VisMontage.h;

% what action
switch (varargin{3}(:)')

    % handle checkbox
    case {'chbox'}

        % if checked
        if tags.CB_vismontage_anatransp.Value > 0

            % set color button enabled
            tags.BT_vismontage_anabackgc.Enable = 'on';

        % otherwise
        else

            % disabled
            tags.BT_vismontage_anabackgc.Enable = 'off';
        end

    % handle color picker
    case {'pick'}

        % show colorpicker
        try
            tags.BT_vismontage_anabackgc.CData = repmat(shiftdim(uint8( ...
                colorpicker(double(lsqueeze( ...
                tags.BT_vismontage_anabackgc.CData(1, 1, :))'), ...
                'background')), -1), 12, 36);
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end
end
