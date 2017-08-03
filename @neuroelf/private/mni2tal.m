function talout = mni2tal(mniin, rounded)
% mni2tal  - converts coordinates from MNI brain to best Talairach guess
%
% FORMAT:       talout = mni2tal(mniin [, rounded])
%
% Input fields:
%
%       mniin       N-by-3 or 3-by-N matrix of coordinates
%       rounded     1x1 double, if given, coordinates are rounded to
%                   specified number of digits
%
% Output fields:
%
%       talout      is the coordinate matrix with Talairach points
%
% (c) Matthew Brett 10/8/99
%
% See also spmtrf

% Version:  v0.9a
% Build:    10051716
% Date:     May-17 2010, 10:48 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% argument check
if nargin < 1 || ...
    length(size(mniin)) > 2 || ...
   ~isa(mniin, 'double')
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing input argument mniin.' ...
    );
end

dimdim = find(size(mniin) == 3);
if isempty(dimdim)
    error( ...
        'neuroelf:BadArguments', ...
        'mniin argument must be a N-by-3 or 3-by-N matrix' ...
    );
end

% transpose as needed
if dimdim(1) == 2
    mniin = mniin';
end

% Transformation matrices, different zooms above/below AC
upT   = spmtrf([0, 0, 0], [0.05, 0, 0], [0.99, 0.97, 0.92]);
downT = spmtrf([0, 0, 0], [0.05, 0, 0], [0.99, 0.97, 0.84]);

% find points below axial AC plane
tmp = mniin(3,:) < 0;  % 1 if below AC

% add 1's for matrix multiplication
talout = [mniin; ones(1, size(mniin, 2))];

% multiply according to above/below axial AC plane
talout(:,  tmp) = downT * talout(:,  tmp);
talout(:, ~tmp) = upT   * talout(:, ~tmp);

% return only 3 coordinate points
talout = talout(1:3, :);

% retranspose to match input
if dimdim(1) == 2
    talout = talout';
end

% round output
if nargin > 1 && ...
    isa(rounded, 'double') && ...
   ~isempty(rounded)
    talout = (1/10^rounded(1)) * round(talout * 10^rounded(1));
end
