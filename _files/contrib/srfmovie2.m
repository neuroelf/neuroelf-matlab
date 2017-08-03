% clear everything
fclose all;
xfig = xfigure;
xfig.DeleteAllFigures;
clear classes;

% change into folder
cd /Data/guthead

% find subject and group folders
gths = findfiles(pwd, 'gth*', 'dirs', 'depth=1', 'relative=');

% find images
pngs = cell(18, 4);
for c = 1:18
    pngs{c, 1} = findfiles([pwd '/' gths{c}], '*LH*lat*.png');
    pngs{c, 2} = findfiles([pwd '/' gths{c}], '*LH*med*.png');
    pngs{c, 3} = findfiles([pwd '/' gths{c}], '*RH*med*.png');
    pngs{c, 4} = findfiles([pwd '/' gths{c}], '*RH*lat*.png');
end

% create image
im = uint8(zeros(720, 1280, 3));
imm = uint8(zeros(720, 256, 3));

% iterate over time (15 seconds by 24 frames + 1)
ip = cell(18, 4);
for tc = 1:361

    % read in the required images
    for ic = 1:72
        ip{ic} = imread(pngs{ic}{tc});
    end

    % stack images
    ipl1 = cat(1, ip{1:8, 1});
    ipl2 = cat(1, ip{1:8, 2});
    ipl3 = cat(1, ip{1:8, 3});
    ipl4 = cat(1, ip{1:8, 4});
    ipr1 = cat(1, ip{10:17, 1});
    ipr2 = cat(1, ip{10:17, 2});
    ipr3 = cat(1, ip{10:17, 3});
    ipr4 = cat(1, ip{10:17, 4});

    % set left and right parts
    im(9:712, 1:512, :) = cat(2, ipl1, ipl2, ipl3, ipl4);
    im(9:712, 513:768, :) = cat(1, ip{18, :});
    im(9:712, 769:1280, :) = cat(2, ipr1, ipr2, ipr3, ipr4);

    % store
    imwrite(im, sprintf('%s/cgth/cgth_frame%03d.png', pwd, tc));
end
