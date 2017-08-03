% example file using SVM on average-extracts from VOI betas (2nd level)

% import from neuroelf library
using(neuroelf, {'bvcoordconv', 'histcount', 'ne_svmpredict', 'ne_svmtrain', 'scaleimage'});

% settings
% -> percentage of cases used for training
trainpct = .75;

% -> iterations (outer and inner)
niter = 200;
lniter = 200;

% -> which maps are used for feature extraction
maps = [3];

% indices of controls and patients
controls = 1:22;
patients = 23:44;

% generate labels
labels = [ones(numel(controls), 1); 2 .* ones(numel(patients), 1)];

% and selection
csubsel = floor(trainpct * numel(controls));
psubsel = floor(trainpct * numel(patients));

% access GLM (re-write to load specific GLM!)
x = xff;
glms = x.Documents('glm');
if numel(glms) ~= 1
    error('more than one GLM loaded.');
end
glm = x.Document(glms{1});
no_subjects = numel(glm.Subjects);
if no_subjects ~= (numel(controls) + numel(patients))
    error('number of subjects mismatch.');
end

% access VOI
vois = x.Documents('voi');
if numel(vois) ~= 1
    error('more than one VOI loaded.');
end
voi = x.Document(vois{1});

% access the VOI voxels
voivox = voi.VOI;
voivox = cat(1, voivox.Voxels);
voigidx = bvcoordconv(voivox, 'tal2bvx', glm.BoundingBox);

% compute histogram
voihist = histcount(voigidx, 1, max(voigidx)+1, 1);

% only select voxels that are at least 25% covered
voigidx = find(voihist > 0);

% instantiate feature space
no_features = numel(maps) * numel(voigidx);
features = zeros(no_subjects, no_features);

% fill features from GLM data
for sc = 1:no_subjects

    % access maps
    for mc = 1:numel(maps)
        map = glm.GLMData.Subject(sc).BetaMaps(:, :, :, maps(mc));
        features(sc, ((mc-1) * numel(voigidx) + 1):(mc * numel(voigidx))) = map(voigidx);
    end
end

% show features (check!)
% scaleimage(features);
drawnow;

% iterate 1000 times
class_accuracy = zeros(niter, lniter);
class_labels = zeros(numel(labels) - (csubsel + psubsel), niter, lniter);
for lc = 1:lniter
    for ic = 1:niter

        % draw up to 80% of data for training
        cselect = randperm(numel(controls));
        pselect = numel(controls) + randperm(numel(patients));
        ctrain = cselect(1:csubsel);
        cclass = cselect(csubsel+1:end);
        ptrain = pselect(1:psubsel);
        pclass = pselect(psubsel+1:end);

        % combine
        trainidx = [ctrain(:)', ptrain(:)'];
        classidx = [cclass(:)', pclass(:)'];

        % train model
        model = ne_svmtrain(labels(trainidx), features(trainidx, :));

        % test model
        [predlabs, predacc] = ne_svmpredict(labels(classidx), features(classidx, :), model);

        % store result
        class_accuracy(ic, lc) = predacc(1);
        class_labels(:, ic, lc) = predlabs(:);
    end

    % reshuffle labels
    labels = labels(randperm(numel(labels)));
    disp(lc); pause(0.01);
end
