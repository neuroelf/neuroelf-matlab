clear; clc;

%insert into fmr & stc
NrOfRows = 64;
NrOfCols = 64;
NVols = 50;
NSlices = 10;


%% create new template fmr
fmr = xff('new:fmr');
fmr.NrOfVolumes = NVols;
fmr.NrOfSlices = NSlices;
fmr.Slice.STCData = uint16(zeros(NrOfRows, NrOfCols, NVols, NSlices));
fmr.SaveAs('Test.fmr');

% fmr sanity check
fmr
fmr.clearobject;
fmr2 = xff('Test.fmr')
fmr2

%% create new stc template
stc = xff('new:stc');
stc.NrOfRows = NrOfRows;
stc.NrOfCols = NrOfCols;
stc.NrOfVolumes = NVols;
stc.NrOfSlices = NSlices;
stc.STCData = zeros(NrOfRows, NrOfCols, NVols, NSlices);
stc.SaveAs('Test.stc');

% stc sanity check
stc
stc.clearobject;
stc2 = xff('Test.stc');
stc2

fmr2.clearobject;
stc2.clearobject;



