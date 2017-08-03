% clear everything
fclose all;
xfig = xfigure;
xfig.DeleteAllFigures;
clear classes;

% change directory
cd /Data/guthead

% open GUI
neuroelf_gui;

% get handles
global ne_gcfg;
ch = ne_gcfg.h;

% load surfaces
lh = xff('LH_smoothed.srf');
rh = xff('RH_smoothed.srf');

% get surface files
lms = findfiles(pwd, '*LH.mtc');
rms = findfiles(pwd, '*RH.mtc');

% get subject IDs
gths = findfiles(pwd, 'gth*', 'dirs', 'depth=1', 'relative=');

% browse left surface
lh.Browse;

% iterate over subjects (and group)
for sc = 18:18

    % load LH MTC
    m = xff(lms{sc});

    % settings
    m.Map(1).ShowPositiveNegativeFlag = 1;
    m.Map(2).ShowPositiveNegativeFlag = 1;
    if sc < 18
        m.ConditionThresholds(:, 2, 1) = 2.5;
        m.ConditionThresholds(:, 2, 2) = 5;
    else
        m.ConditionThresholds(:, 2, 1) = 4.5;
        m.ConditionThresholds(:, 2, 2) = 9;
    end

    % iterate over time (one volume is 500, so we use 1/12 timing!)
    for tc = 1:(1/12):41

        % browse
        m.Browse(1:2, tc);

        % set surface position
        neuroelf_gui('setsurfpos', [], {180, 0, [0, -25 -20], 4/3});

        % make sure figure is visible
        figure(ch.MainFigMLH);
        drawnow;

        % get frame from axes
        fr = getframe(ch.Surface);

        % resize portion
        if sc < 18
            im = image_resize(fr.cdata(77:428, 1:512, :), 88, 128);
        else
            im = image_resize(fr.cdata(77:428, 1:512, :), 176, 256);
        end

        % save as frame
        imwrite(im, sprintf('%s/%s_LH_lat_%05dms.png', gths{sc}, gths{sc}, ...
            round((tc-1) * 500)));

        % set other side
        neuroelf_gui('setsurfpos', [], {0, 0, [0, 24 -20], 4/3});

        % make sure figure is visible
        figure(ch.MainFigMLH);
        drawnow;

        % get frame from axes
        fr = getframe(ch.Surface);

        % resize portion
        if sc < 18
            im = image_resize(fr.cdata(77:428, 1:512, :), 88, 128);
        else
            im = image_resize(fr.cdata(77:428, 1:512, :), 176, 256);
        end

        % save as frame
        imwrite(im, sprintf('%s/%s_LH_med_%05dms.png', gths{sc}, gths{sc}, ...
            round((tc-1) * 500)));
    end

    % clear MTC
    m.ClearObject;
end

% unbrowse left and browse right surface instead
lh.ClearObject;
rh.Browse;

% iterate over subjects (and group)
for sc = 1:18

    % load RH MTC
    m = xff(rms{sc});

    % settings
    m.Map(1).ShowPositiveNegativeFlag = 1;
    m.Map(2).ShowPositiveNegativeFlag = 1;
    if sc < 18
        m.ConditionThresholds(:, 2, 1) = 2.5;
        m.ConditionThresholds(:, 2, 2) = 5;
    else
        m.ConditionThresholds(:, 2, 1) = 4.5;
        m.ConditionThresholds(:, 2, 2) = 9;
    end

    % iterate over time (one volume is 500, so we use 1/12 timing!)
    for tc = 1:(1/12):33

        % browse
        m.Browse(1:2, tc);

        % set surface position
        neuroelf_gui('setsurfpos', [], {180, 0, [0, -25 -20], 4/3});

        % make sure figure is visible
        figure(ch.MainFigMLH);
        drawnow;

        % get frame from axes
        fr = getframe(ch.Surface);

        % resize portion
        if sc < 18
            im = image_resize(fr.cdata(77:428, 1:512, :), 88, 128);
        else
            im = image_resize(fr.cdata(77:428, 1:512, :), 176, 256);
        end

        % save as frame
        imwrite(im, sprintf('%s/%s_RH_med_%05dms.png', gths{sc}, gths{sc}, ...
            round((tc-1) * 500)));

        % set other side
        neuroelf_gui('setsurfpos', [], {0, 0, [0, 24 -20], 4/3});

        % make sure figure is visible
        figure(ch.MainFigMLH);
        drawnow;

        % get frame from axes
        fr = getframe(ch.Surface);

        % resize portion
        if sc < 18
            im = image_resize(fr.cdata(77:428, 1:512, :), 88, 128);
        else
            im = image_resize(fr.cdata(77:428, 1:512, :), 176, 256);
        end

        % save as frame
        imwrite(im, sprintf('%s/%s_RH_lat_%05dms.png', gths{sc}, gths{sc}, ...
            round((tc-1) * 500)));
    end

    % clear MTC
    m.ClearObject;
end

% clear surface
rh.ClearObject;
