function multiset(handles,varargin)
% multiset  - set properties for multiple handles
%
% FORMAT:       multiset(handles, property, values, ...)
%       OR      multiset(handles, propstruct)
%
% Input fields:
%
%       handles     Nx1 double array with object handles
%       property    string naming the property to set
%       values      Nx1 cell array with values to set
%       propstruct  Nx1 struct array with properties/values
%
% See also get, subget, set.

% Version:  v0.9b
% Build:    11050712
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

% argument check
if nargin < 2
    error( ...
        'neuroelf:TooFewArguments', ...
        'Too few arguments. Try ''help %s''.', ...
        mfilename ...
    );
end
if nargin == 2 && ...
   ~isstruct(varargin{1})
    error( ...
        'neuroelf:BadArgument', ...
        'No properties/values specified.' ...
    );
end
if ~isa(handles, 'double')
    error( ...
        'neuroelf:BadArgument', ...
        'Wrong input type for handles.' ...
    );
end

% return on no handles
if isempty(handles)
    return;
end

% get number of handles
nh = numel(handles);

% check further for struct
if isstruct(varargin{1})
    ps = varargin{1};
    ns = numel(ps);
    if ns < nh && ...
        ns ~= 1
        error( ...
            'neuroelf:BadArgument', ...
            'Not enough property values specified.' ...
        );
    end
    if ns == 1 && ...
        nh > 1
        ps(1:nh) = ps;
    end
    pf = fieldnames(ps);
    nf = length(pf);
    for fc = 1:nf
        try
            multiset(handles,pf{fc},{ps(:).(pf{fc})}');
        catch ne_eo;
            warning( ...
                'neuroelf:InvalidProperty', ...
                'Couldn''t set property ''%s'' (%s)', ...
                pf{fc}, ne_eo.message ...
            );
        end
    end
else
    for vac = 1:2:(nargin-1)
        if ~isrealvarname(varargin{vac})
            error( ...
                'neuroelf:InvalidProperty', ...
                'Invalid property name given.' ...
            );
        end
        if ~iscell(varargin{vac+1})
            error( ...
                'neuroelf:InvalidPropertyValue', ...
                'Invalid property value(s) given.' ...
            );
        end
        pn = varargin{vac};
        pv = varargin{vac+1};
        nv = numel(pv);
        if nv < nh && ...
            nv ~= 1
            error( ...
                'neuroelf:BadArgument', ...
                'Not enough property values specified.' ...
            );
        end
        if nv == 1 && ...
            nh > 1
            pv = pv(ones(1,nh));
        end
        for hc = 1:nh
            try
                set(handles(hc),pn,pv{hc});
            catch ne_eo;
                error( ...
                    'neuroelf:InvalidProperty', ...
                    'Error setting property. (%s)', ...
                    ne_eo.message ...
                );
            end
        end
    end
end
