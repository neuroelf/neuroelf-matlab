% import from neuroelf library
using(neuroelf, {'alphasim', 'findfirst', 'lsqueeze', 'newnatresvmp'});

% load mask file (66-by-52-by-56)
mask = xff([neuroelf_path('masks') '/colin_brain_ICBMnorm_wholebrain3mm.msk']);
m = (mask.Mask > 0);

% number of iterations, number of maps to simulate
niter = 10000;
tdf = 30;

% list of thresholds
thr = [0.01 .* (5:-0.5:2), 0.001 .* (19:-1:1), 0.0001 .* (9:-1:1), 0.00001 .* (9:-1:1), 1e-6 .* (9:-1:1), 1e-7 .* (9:-2:1), 1e-8 .* (9:-2:1), 1e-9 .* [5, 2, 1]];

% list of z-shifts
zshifts = lsqueeze([0:0.05:1.25; 0:-.05:-1.25]);

% list of FWHM values
fwhm = [1:0.05:2, 2.1:0.1:4, 4.2:.2:7];

% tails
tails = {[2, 2], [1, 2]};

% generate VMP?
vmpfile = sprintf('%s/alphasim_table_%diters_%dtdf.vmp', neuroelf_path('contrib'), niter, tdf);
if exist(vmpfile, 'file') ~= 2

    % new VMP
    vmp = newnatresvmp();

    % make settings
    vmp.Resolution = 3;
    vmp.Map.VMPData = single(zeros(size(m)));
    vmp.Map.Type = 55;
    vmp.Map.Name = sprintf('alphasim table (%d of %d tails)', tails{1}(1), tails{1}(2));
    vmp.Map.LowerThreshold = 2;
    vmp.Map.UpperThreshold = 100;
    vmp.Map = vmp.Map([1, 1]);
    vmp.Map(2).Name = sprintf('alphasim table (%d of %d tails)', tails{2}(1), tails{2}(2));
    %vmp.Map(3).Name = sprintf('alphasim table (%d of %d tails)', tails{3}(1), tails{3}(2));
    vmp.NrOfMaps = 3;
    vmp.XStart = mask.XStart;
    vmp.XEnd = mask.XEnd;
    vmp.YStart = mask.YStart;
    vmp.YEnd = mask.YEnd;
    vmp.ZStart = mask.ZStart;
    vmp.ZEnd = mask.ZEnd;

    % additional settings
    vmp.RunTimeVars.AutoSave = true;
    vmp.RunTimeVars.FWHM = fwhm;
    vmp.RunTimeVars.Mask = uint8(m);
    vmp.RunTimeVars.NIter = niter;
    vmp.RunTimeVars.Thresholds = thr;
    vmp.RunTimeVars.StatsType = tails;
    vmp.RunTimeVars.zShift = zshifts;

    % save VMP
    vmp.SaveAs(vmpfile);

    % clear
    vmp.ClearObject;
end

% load VMP
vmp = bless(xff(vmpfile, 'T'));

% get data from VMP
fwhm = vmp.RunTimeVars.FWHM;
niter = vmp.RunTimeVars.NIter;
tails = vmp.RunTimeVars.StatsType;
thr = vmp.RunTimeVars.Thresholds;
zshifts = vmp.RunTimeVars.zShift;

% iterate over z-shifts
for zc = 2:52

    % iterate over FWHM values
    for fc = 1:56

        % iterate over tails/stype
        for tc = 1:numel(tails)

            % skip if done
            if any(vmp.Map(tc).VMPData(:, zc, fc) ~= 0)
                continue;
            end

            % mark as being done
            vmp.Map(tc).VMPData(66, zc, fc) = -1;
            disp([tc, NaN, zc, fc, double(any(vmp.Map(tc).VMPData(:, zc, fc) ~= 0))]);
            drawnow;
            pause(0.01);

            % run alphasim
            a = alphasim(size(m), struct( ...
                'fwhm',  fwhm(fc) .* ones(1, 3), ...
                'mask',   m, ...
                'niter',  niter, ...
                'tdf',    tdf, ...
                'thr',    thr, ...
                'stype',  tails{tc}, ...
                'zshift', zshifts(zc)));

            % fill in data
            for thc = 1:66

                % find required minimal cluster size
                mcs = findfirst(a{thc}(:, end) >= 0.05, -1);

                % store
                if zc == 2
                    vmp.Map(tc).VMPData(thc, 1, fc) = mcs;
                end
                vmp.Map(tc).VMPData(thc, zc, fc) = mcs;
            end
        end
    end
end
