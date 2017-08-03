function mniout = tal2mni(talin, rounded)
% tal2mni  - converts coordinates from TAL brain to best MNI-space guess
%
% FORMAT:       mniout = tal2mni(talin [, rounded])
%
% Input fields:
%
%       talin       N-by-3 or 3-by-N matrix of coordinates
%       rounded     1x1 double, if given, coordinates are rounded to
%                   specified number of digits
%
% Output fields:
%
%       mniout      is the coordinate matrix with MNI-space points
%
% (c) Matthew Brett 2/2/01
%
% See also spmtrf

% Version:  v0.9a
% Build:    10051716
% Date:     May-17 2010, 10:48 AM EST
% Editor:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% persistent transformation matrices
persistent t2m_trf;
if isempty(t2m_trf) || ...
   ~isstruct(t2m_trf)
    t2m_trf = struct( ...
        'rotT',  inv(spmtrf([0, 0, 0], [0.05, 0, 0], [1, 1, 1])), ...
        'upT',   inv(spmtrf([0, 0, 0], [0, 0, 0], [0.99, 0.97, 0.92])), ...
        'downT', inv(spmtrf([0, 0, 0], [0, 0, 0], [0.99, 0.97, 0.84])));
end

% argument check
if nargin < 1 || ...
    length(size(talin)) > 2 || ...
   ~isa(talin, 'double')
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing input argument mniin.' ...
    );
end

% transpose ?
dimdim = find(size(talin) == 3);
if isempty(dimdim)
    error( ...
        'neuroelf:BadArguments', ...
        'talin argument must be a N-by-3 or 3-by-N matrix' ...
    );
end

% transpose as needed
if dimdim(1) == 2
    talin = talin';
end

% add row of ones for multiplication
talin = [talin; ones(1, size(talin, 2))];

% "un-rotate" first
talin = t2m_trf.rotT * talin;

% find indices below Z == 0
tmp = talin(3, :) < 0;

% multiply according to above/below axial AC plane
talin(:,  tmp) = t2m_trf.downT * talin(:,  tmp);
talin(:, ~tmp) = t2m_trf.upT   * talin(:, ~tmp);

% return only 3 coordinate points
mniout = talin(1:3, :);

% retranspose to match input
if dimdim(1) == 2
    mniout = mniout';
end

% round output
if nargin > 1 && ...
    isa(rounded, 'double') && ...
   ~isempty(rounded)
    mniout = (1/10^rounded(1)) * round(mniout * 10^rounded(1));
end
