function [label, center, bCon, sumD, D] = litekmeans(X, k, varargin)
% litekmeans  - K-means clustering, accelerated by matlab matrix operations
%
% FORMAT:       [l, c, bc, sd, d] = litekmeans(x, k, ['PARAM', value, ...])
%
% Input fields:
%
%       x           SxD double matrix (samples-by-dimensions)
%       k           number of class labels to produce
%                   supported parameter/value pairs are
%       'ClusterMaxIter'  - number of iterations for seed finding (10)
%       'Distance'  distance measure (type) that will be minimzed, one of
%         {'sqEuclidean'} - squared Euclidean distance
%          'cosine'       - 1 minus the cosine of the included angle
%                           between points (treated as vectors); each row
%                           of x should be normalized to unit; if the
%                           intial center matrix is provided, it should
%                           also be normalized
%       'MaxIter'   maximum number of iterations allowed (default: 100)
%       'Replicates number of times to repeat the clustering, each with a
%                   new set of initial centroids (default: 1; if the initial
%                   centroids are provided, the replicate will be set to 1)
%       'Start'     method used to choose initial cluster centroids
%         {'sample'} - Select K observations from X at random
%          'cluster' - Perform preliminary clustering phase on random 10%
%                      subsample of X; this preliminary phase is itself
%                      initialized using 'sample'; An additional parameter,
%                      'clusterMaxIter', can be used to control the maximum
%                      number of iterations in each preliminary clustering
%                      problem
%           matrix   - A K-by-P matrix of starting locations; or a K-by-1
%                      indicate vector indicating which K points in X
%                      should be used as the initial center; in this case,
%                      you can pass in [] for K, and KMEANS infers K from
%                      the first dimension of the matrix
%
% Output fields:
%
%       l           Sx1 class/cluster labels (range 1:k)
%       c           KxD class/cluster centroid points
%       bc          boolean value indicating whether the iteration converged
%       sd          1xK within-cluster sums of distances
%       d           SxK distances of each point to each clusters centroid
%
% Note: the algorithm partitions the points in the S-by-D data matrix, x,
%       into k clusters; this partition minimizes the sum, over all
%       clusters, of the within-cluster sums of point-to-cluster-centroid
%       distances
%       rows of x correspond to points, columns correspond to variables
%       litekmeans returns an S-by-1 vector label containing the
%       cluster indices of each sample/point
%
% Examples:
%
%       fea = rand(500,10);
%       [label, center] = litekmeans(fea, 5, 'MaxIter', 50);
%
%       fea = rand(500,10);
%       [label, center] = litekmeans(fea, 5, 'MaxIter', 50, 'Replicates', 10);
%
%       fea = rand(500,10);
%       [label, center, bCon, sumD, D] = litekmeans(fea, 5, 'MaxIter', 50);
%       TSD = sum(sumD);
%
%       fea = rand(500,10);
%       initcenter = rand(5,10);
%       [label, center] = litekmeans(fea, 5, 'MaxIter', 50, 'Start', initcenter);
%
%       fea = rand(500,10);
%       idx=randperm(500);
%       [label, center] = litekmeans(fea, 5, 'MaxIter', 50, 'Start', idx(1:5));
%
% See also KMEANS
%
% [Cite] Deng Cai, "Litekmeans: the fastest matlab implementation of
%        kmeans" - available at:
%        http://www.zjucadcg.cn/dengcai/Data/Clustering.html, 2011.

% Version:  v2.0
% Build:    13020213
% Date:     Feb-02 2013, 1:49 PM EST
% Author:   Deng Cai (dengcai AT gmail.com)
% Editor:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, Deng Cai
% All rights reserved.

% version 2.0 --December/2011
% version 1.0 --November/2011

% argument check
if nargin < 2 || ...
   ~isa(X, 'double') || ...
    ndims(X) ~= 2 || ...
    numel(k) ~= 1
    error( ...
        'litekmeans:TooFewInputs', ...
        'At least two input arguments, x and k, required.' ...
    );
end
[n, p] = size(X);

% parse optional arguments
pnames = {   'distance', 'start' ,  'maxiter',  'replicates', 'onlinephase', 'clustermaxiter'};
dflts =  {'sqeuclidean', 'sample',      []   ,       []     ,         'off',        []       };
[eid, errmsg, distance, start, maxit, reps, online, clustermaxit] = ...
    getargs(pnames, dflts, varargin{:});
if ~isempty(eid)
    error(sprintf('litekmeans:%s', eid), errmsg);
end
if ischar(distance)
    distNames = {'sqeuclidean','cosine'};
    j = strcmpi(distance, distNames);
    j = find(j);
    if length(j) > 1
        error('litekmeans:AmbiguousDistance', ...
            'Ambiguous ''Distance'' parameter value:  %s.', distance);
    elseif isempty(j)
        error('litekmeans:UnknownDistance', ...
            'Unknown ''Distance'' parameter value:  %s.', distance);
    end
    distance = distNames{j};
else
    error('litekmeans:InvalidDistance', ...
        'The ''Distance'' parameter value must be a string.');
end
center = [];
if ischar(start)
    startNames = {'sample', 'cluster'};
    j = find(strncmpi(start, startNames, length(start)));
    if length(j) > 1
        error(message('litekmeans:AmbiguousStart', start));
    elseif isempty(j)
        error(message('litekmeans:UnknownStart', start));
    elseif isempty(k)
        error('litekmeans:MissingK', ...
            'You must specify the number of clusters, K.');
    end
    if j == 2
        if floor(.1 * n) < 5 * k
            j = 1;
        end
    end
    start = startNames{j};
elseif isnumeric(start)
    if size(start, 2) == p
        center = start;
    elseif (size(start, 2) == 1 || size(start, 1) == 1)
        center = X(start, :);
    else
        error('litekmeans:MisshapedStart', ...
            'The ''Start'' matrix must have the same number of columns as X.');
    end
    if isempty(k)
        k = size(center, 1);
    elseif (k ~= size(center, 1))
        error('litekmeans:MisshapedStart', ...
            'The ''Start'' matrix must have K rows.');
    end
    start = 'numeric';
else
    error('litekmeans:InvalidStart', ...
        'The ''Start'' parameter value must be a string or a numeric matrix or array.');
end

% the maximum iteration number is default 100
if isempty(maxit)
    maxit = 100;
end

% the maximum iteration number for preliminary clustering phase on random
% 10% subsamples is default 10
if isempty(clustermaxit)
    clustermaxit = 10;
end

% assume one replicate
if isempty(reps) || ...
   ~isempty(center)
    reps = 1;
end
if ~(isscalar(k) && isnumeric(k) && isreal(k) && k > 0 && (round(k) == k))
    error('litekmeans:InvalidK', ...
        'X must be a positive integer value.');
elseif n < k
    error('litekmeans:TooManyClusters', ...
        'X must have more rows than the number of clusters.');
end

% initialize some outputs
bestlabel = [];
sumD = zeros(1, k);
bCon = false;

% iterate across replications
for t = 1:reps

    % permute points
    rs = randperm(n);

    % centroid selection
    switch start

        % from data
        case 'sample'
            center = X(rs(1:k), :);

        % cluster instead
        case 'cluster'
            [dump, center] = litekmeans(X(rs(ceil(.1 * n)), :), k, varargin{:}, ...
                'start', 'sample', 'replicates', 1, 'MaxIter', clustermaxit);

        % data given already!
        case 'numeric'
    end

    % initialize counters
    last = 0;
    label = 1;
    it = 0;

    % distance measure
    switch distance

        % euclidean
        case 'sqeuclidean'

            % find labels (up to max. iterations)
            while any(label ~= last) && it < maxit

                % keep track of labels
                last = label;

                % distance computation
                bb = full(sum(center .* center, 2)');
                ab = full(X * center');
                D = bb(ones(1, n), :) - 2 * ab;

                % assign samples to the nearest centers
                [val, label] = min(D, [], 2);
                ll = unique(label);

                % some clusters dropped ...
                if length(ll) < k
                    missCluster = 1:k;
                    missCluster(ll) = [];
                    missNum = length(missCluster);
                    aa = sum(X .* X, 2);
                    val = aa + val;
                    [dump, idx] = sort(val, 1, 'descend');
                    label(idx(1:missNum)) = missCluster;
                end

                % transform label into indicator matrix
                E = sparse(1:n,label,1,n,k,n);

                % compute center of each cluster
                center = full((E * spdiags(1 ./ sum(E, 1)', 0, k, k))' * X);

                % iteration counter
                it = it + 1;
            end

            % converged?
            if it < maxit
                bCon = true;
            end

            % assign best labels
            if isempty(bestlabel)
                bestlabel = label;
                bestcenter = center;

                % with replications?
                if reps > 1

                    % not converged
                    if it >= maxit
                        aa = full(sum(X.*X,2));
                        bb = full(sum(center.*center,2));
                        ab = full(X*center');
                        D = bsxfun(@plus,aa,bb') - 2*ab;
                        D(D<0) = 0;

                    % converged
                    else
                        aa = full(sum(X.*X,2));
                        D = aa(:,ones(1,k)) + D;
                        D(D<0) = 0;
                    end

                    % distance measure
                    D = sqrt(D);

                    % compute across clusters
                    for j = 1:k
                        sumD(j) = sum(D(label==j,j));
                    end
                    bestsumD = sumD;
                    bestD = D;
                end

            % update best labels
            else

                % not converged
                if it >= maxit
                    aa = full(sum(X.*X,2));
                    bb = full(sum(center.*center,2));
                    ab = full(X*center');
                    D = bsxfun(@plus,aa,bb') - 2*ab;
                    D(D<0) = 0;

                % converged
                else
                    aa = full(sum(X.*X,2));
                    D = aa(:,ones(1,k)) + D;
                    D(D<0) = 0;
                end

                % distance measure
                D = sqrt(D);
                for j = 1:k
                    sumD(j) = sum(D(label==j,j));
                end
                if sum(sumD) < sum(bestsumD)
                    bestlabel = label;
                    bestcenter = center;
                    bestsumD = sumD;
                    bestD = D;
                end
            end

        % cosine distance
        case 'cosine'

            % follows logic from above
            while any(label ~= last) && it < maxit
                last = label;

                % compute measure and assign samples to nearest center
                W = full(X * center');
                [val, label] = max(W, [], 2);
                ll = unique(label);
                if length(ll) < k
                    missCluster = 1:k;
                    missCluster(ll) = [];
                    missNum = length(missCluster);
                    [dump, idx] = sort(val);
                    label(idx(1:missNum)) = missCluster;
                end

                % transform label into indicator matrix
                E = sparse(1:n, label, 1, n, k, n);

                % compute center of each cluster
                center = full((E * spdiags(1 ./ sum(E, 1)', 0, k, k))' * X);
                centernorm = sqrt(sum(center.^2, 2));
                center = center ./ centernorm(:, ones(1, p));

                % iteration counter
                it = it + 1;
            end
            if it<maxit
                bCon = true;
            end
            if isempty(bestlabel)
                bestlabel = label;
                bestcenter = center;
                if reps > 1
                    if any(label ~= last)
                        W = full(X * center');
                    end
                    D = 1 - W;
                    for j = 1:k
                        sumD(j) = sum(D(label == j, j));
                    end
                    bestsumD = sumD;
                    bestD = D;
                end
            else
                if any(label ~= last)
                    W=full(X * center');
                end
                D = 1 - W;
                for j = 1:k
                    sumD(j) = sum(D(label == j, j));
                end
                if sum(sumD) < sum(bestsumD)
                    bestlabel = label;
                    bestcenter = center;
                    bestsumD = sumD;
                    bestD = D;
                end
            end
    end
end

% fill outputs
label = bestlabel;
center = bestcenter;
if reps > 1
    sumD = bestsumD;
    D = bestD;
elseif nargout > 3
    switch distance

        % euclidean distance measure
        case 'sqeuclidean'
            aa = full(sum(X .* X, 2));
            if it >= maxit
                bb = full(sum(center .* center, 2));
                ab = full(X * center');
                D = bsxfun(@plus, aa, bb') - 2 * ab;
            else
                D = aa(:, ones(1, k)) + D;
            end
            D(D < 0) = 0;
            D = sqrt(D);

        % cosine distance measure
        case 'cosine'
            if it >= maxit
                W = full(X * center');
            end
            D = 1 - W;
    end
    for j = 1:k
        sumD(j) = sum(D(label == j, j));
    end
end


% sub-function for parsing arguments (modified from MATLAB)
function [eid, emsg, varargout] = getargs(pnames, dflts, varargin)
% original Copyright 1993-2008 The MathWorks, Inc.
% modified by Deng Cai (dengcai@gmail.com) 2011.11.27

% initialize some variables
emsg = '';
eid = '';
nparams = numel(pnames);
varargout = dflts;
unrecog = {};
nargs = numel(varargin);

% must have name/value pairs
if mod(nargs, 2)~=0
    eid = 'WrongNumberArgs';
    emsg = 'Wrong number of arguments.';
else

    % process name/value pairs
    for j = 1:2:nargs
        pname = varargin{j};
        if ~ischar(pname)
            eid = 'BadParamName';
            emsg = 'Parameter name must be text.';
            break;
        end
        i = strcmpi(pname, pnames);
        i = find(i);
        if isempty(i)
            % if they've asked to get back unrecognized names/values, add this
            % one to the list
            if nargout > nparams+2
                unrecog((end+1):(end+2)) = {varargin{j}, varargin{j+1}};
                % otherwise, it's an error
            else
                eid = 'BadParamName';
                emsg = sprintf('Invalid parameter name:  %s.', pname);
                break;
            end
        elseif length(i) > 1
            eid = 'BadParamName';
            emsg = sprintf('Ambiguous parameter name:  %s.', pname);
            break;
        else
            varargout{i} = varargin{j+1};
        end
    end
end

% return everything
varargout{nparams+1} = unrecog;
