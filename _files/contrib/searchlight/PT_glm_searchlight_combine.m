function varargout = PT_glm_searchlight_combine(slfolder)
%PT_glm_searchlight_combine  Combine sub-VMPs into a single one.
%   This function does not return any outputs. If no input is given, it
%   will ask the user to specify a folder name.
%
%   The target VMP will be called after the last fileparts part of the
%   input folder + .vmp. All partial VMPs (*_vmp.mat) will be deleted after
%   they have been successfully processed.

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% input
if nargin < 1 || ~ischar(slfolder) || exist(slfolder(:)', 'dir') ~= 7
    slfolder = uigetdir(pwd, 'Please choose a folder containing the SL data...');
    if isequal(slfolder, 0) || isempty(slfolder) || ~ischar(slfolder) || exist(slfolder, 'dir') ~= 7
        return;
    end
end
slfolder = slfolder(:)';
[~, slfile] = fileparts(slfolder);
if isempty(slfile)
    [~, slfile] = fileparts(fileparts(slfolder));
end

% require neuroelf
n = neuroelf;

% find *_*of*_vmp.mat files
vmpmats = n.findfiles(slfolder, '*_*of*_vmp.mat', '-d1');
if isempty(vmpmats)
    disp('No more *_vmp.mat files to process.');
    return;
end

% locate combined file
cvmpfile = [slfolder filesep slfile '.vmp'];
if exist(cvmpfile, 'file') ~= 2

    % copy first file
    matc = load(vmpmats{1});
    vmpc = matc.vmpc;

    % create VMP
    cvmp = xff('new:vmp');

    % set content
    setcont(cvmp, vmpc);

    % but scratch NULL
    cvmp.RunTimeVars.AutoSave = true;
    cvmp.RunTimeVars.SLNull = zeros(0, 1);
else
    % load VMP
    cvmp = xff(cvmpfile);
end

% mask
if exist([cvmpfile(1:end-4) '_done.msk'], 'file') ~= 2

    % create mask
    msk = xff('new:msk');
    msk.Mask = uint8(zeros(size(cvmp.Map(1).VMPData)));
    msk.Resolution = cvmp.Resolution;
    msk.XStart = cvmp.XStart;
    msk.XEnd = cvmp.XEnd;
    msk.YStart = cvmp.YStart;
    msk.YEnd = cvmp.YEnd;
    msk.ZStart = cvmp.ZStart;
    msk.ZEnd = cvmp.ZEnd;
else
    msk = xff([cvmpfile(1:end-4) '_done.msk']);
end

% iterate over parts
for pc = 1:numel(vmpmats)
    
    % load VMP
    fprintf('Integrating %s...\n', vmpmats{pc});
    matc = load(vmpmats{pc});
    pvmp = matc.vmpc;

    % add NULL data
    cvmp.RunTimeVars.SLNull = [cvmp.RunTimeVars.SLNull; pvmp.RunTimeVars.SLNull];

    % add data
    pmap = pvmp.Map;
    pdata = cat(4, pmap.VMPData);
    psize = size(pdata);
    pdata = reshape(pdata, prod(psize(1:3)), size(pdata, 4));
    pmask = find(sum(pdata ~= 0, 2) > 26);
    plist = find(all(pdata(pmask, :) ~= 0, 1));
    for mc = plist(:)'
        cvmp.Map(mc).VMPData(pmask) = pdata(pmask, mc);
    end
end

% save VMP
cvmp.Save;

% add to "done" mask
pmap = cvmp.Map;
pdata = cat(4, pmap.VMPData);
msk.Mask = uint8(all(pdata(:, :, :, plist) ~= 0, 4));
fprintf('%d voxels done.\n', sum(msk.Mask(:)));
msk.SaveAs([cvmpfile(1:end-4) '_done.msk']);
msk.ClearObject;

% delete combined files
for pc = 1:numel(vmpmats)
    delete(vmpmats{pc});
end
