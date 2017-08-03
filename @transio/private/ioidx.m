function [outs, idx, ridx, cidx, lastcolon] = ioidx(S, vsize)
% ioidx  - get indices to perform IO on
%
% FORMAT:       [outs, idx, ridx, cidx, lastcolon] = ioidx(S, vsize)
%
% Input fields:
%
%       S           1x1 struct as used in subsref / subsasgn
%       vsize       1xN double, variable size to index into
%
% Output fields:
%
%       outs        1xN double, output size
%       idx         indices to address
%       ridx        reverse indices (to finally store into memory)
%       cidx        continuity indexes (within idx!)
%       lastcolon   number of highest index with full indexing range

% Version:  v0.9a
% Build:    10051716
% Date:     May-17 2010, 10:48 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% argument check
if nargin < 2 || ...
   ~isstruct(S) || ...
    numel(S) ~= 1 || ...
   ~isfield(S, 'type') || ...
   ~ischar(S.type) || ...
   ~strcmp(S.type(:)', '()') || ...
   ~isfield(S, 'subs') || ...
   ~isa(vsize, 'double') || ...
    numel(vsize) < 2 || ...
    numel(vsize) > 32 || ...
    numel(vsize) ~= size(vsize, 2) || ...
    any(isinf(vsize) | isnan(vsize) | vsize ~= fix(vsize))
    error( ...
        'transio:BadArgument', ...
        'No or bad S struct or vsize argument given or empty.' ...
    );
end

% get and check S sizes
idx = S.subs(:)';
numidx = numel(idx);
if numidx < 1
    error( ...
        'transio:BadArgument', ...
        'Invalid size requested.' ...
    );
elseif numidx > length(vsize)
    vsize(end+1:numidx) = 1;
end

% initialize reverse indexing array
ridx = cell(1, numidx);
cidx = cell(1, numidx);

% initialize output size, start and end idx
outs = ones(1, numidx);

% keep track of "all-indexing" and continuous subscripts
iscnt = true;
iscol = true;
lastcolon = 0;

% iterate over
for sc = 1:numidx

    % remaining size
    if sc < numidx
        remsize = vsize(sc);
    else
        remsize = prod(vsize(sc:end));
    end

    % chars ...
    if ischar(idx{sc})

        % are only accepted as ':'
        if numel(idx{sc}) ~= 1 || ...
            idx{sc} ~= ':'
            error( ...
                'transio:BadArgument', ...
                'Illegal char indexing in S.' ...
            );
        end

        % as long as prior expressions were complete, update lastcolon
        if iscol
            lastcolon = sc;
        end

        % what size for this index
        outs(sc) = remsize;

        % set ridx, cidx entries
        ridx{sc} = ':';

        % still contiguous?
        if iscnt
            cidx{sc} = [1; remsize + 1];
        else
            cidx{sc} = [];
        end

        % if more than one index, also re-set idx!
        if numidx > 1
            idx{sc} = 1:remsize;
        end

        % don't go on checking
        continue;
    end

    % index size, length and number of elements
    idxsz = size(idx{sc});
    idxln = max(idxsz);
    idxne = prod(idxsz);

    % logical arrays
    if islogical(idx{sc})

        % must be last index, and reject then too large matrices
        if sc < numidx || ...
            idxne > remsize
            error( ...
                'transio:BadArgument', ...
                'Bad positioned or tool large boolean indexing variable.' ...
            );
        end

        % linearize (find creates a 1-D indexing array) and update idxsz
        idx{sc} = find(idx{sc});
        idxsz = length(idx{sc});

    % reject any other non-numerical indexing
    elseif ~isnumeric(idx{sc})
        error( ...
            'transio:BadArgument', ...
            'Only '':'', logical and numerical indexing allowed.' ...
        );

    % reject N-D array if not in last subscript
    elseif sc < numidx && ...
        idxln < idxne && ...
        idxne > 0
        error( ...
            'transio:BadArugment', ...
            'N-D array only allowed as last subscript or too big.' ...
        );
    end

    % get first and last index
    minidx = min(idx{sc}(:));
    maxidx = max(idx{sc}(:));

    % check idx (min, max, argument position, fixed numbers)
    if ~isempty(minidx) && ...
       (minidx < 1 || ...
        maxidx > remsize || ...
        (sc < numidx && maxidx > vsize(sc)) || ...
         any(fix(idx{sc}(:)) ~= idx{sc}(:)))
        error( ...
            'transio:BadArgument', ...
            'Fixed, numerical indexing within bounds required.' ...
        );
    end

    % check for multi dimensional matrices
    if idxln ~= idxne && ...
        idxne > 0

        % set remaining output size to input size
        outs = [outs(1:sc-1), idxsz];

        % and then linearize input
        idx{sc} = idx{sc}(:)';
    end

    % get sorted array and check uniqueness
    [usz, uszi] = unique(idx{sc});
    idxnn = numel(idx{sc});
    uszne = numel(usz);

    % reject double indexing!
    if uszne ~= idxnn

        error( ...
            'transio:BadArgument', ...
            'Index reusing is not supported for IO.' ...
        );

    % full range
    elseif uszne == remsize

        % update lastcolon property
        if iscol
            lastcolon = sc;
        end

        % set output size accordingly
        if idxln == idxne
            outs(sc) = remsize;
        end

        % no sorting required
        if all(usz == idx{sc})
            ridx{sc} = ':';
        else
            ridx{sc} = uszi;
        end

        % contiguity index
        if iscnt
            cidx{sc} = [1; remsize + 1];
        else
            cidx{sc} = [];
        end

    % as long as *something* is requested
    elseif ~isempty(usz)

        % stop full indexing here
        iscol = false;

        % set output size according to requested array
        if idxln == idxne
            outs(sc) = uszne;
        end

        % put into idx
        idx{sc} = usz;

        % reverse indexing needed
        if any(uszi(:) ~= (1:uszne)')
            ridx{sc} = uszi;
        else
            ridx{sc} = ':';
        end

        % contiguous block requested
        if iscnt && ...
           (1 + usz(end) - usz(1)) == uszne

            % take as one block
            cidx{sc} = [1; uszne + 1];

        % prior index contiguous? then find indices where to index reading
        elseif iscnt

            % find contiguous blocks
            cidx{sc} = 1 + [0; find(diff(usz(:)) - 1); uszne];

        % otherwise
        else

            % take each singly
            cidx{sc} = [];
        end

        % stop contiguity tracking
        iscnt = false;

    % nothing requested
    else
        idx{sc} = [];
        cidx{sc} = [];
        ridx{sc} = 1;
    end
end

% needed to get size right!
if numel(outs) < 2
    if isempty(outs)
        outs = [1, 1];
    else
        outs(2) = 1;
    end
end
