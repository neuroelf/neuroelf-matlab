function neuroelf_sptexts(isize)
%NEUROELF_SPTEXTS  Re-create NeuroElf splash image texts.
%   NEUROELF_SPTEXTS(ISIZE) creates the text-based splash images for the
%   specified image size. If no argument is given, it uses the default.

% Version:  v1.1
% Build:    16081613
% Date:     Aug-16 2016, 1:25 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2016, Jochen Weber
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%     * Neither the name of Columbia University nor the
%       names of its contributors may be used to endorse or promote products
%       derived from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% splash image size as default
if nargin < 1 || ~isa(isize, 'double') || numel(isize) ~= 2 || ...
    any(isinf(isize) | isnan(isize) | isize < 1)
    isize = [1248, 448];
end

% generate texts
nv = neuroelf_version;
nb = sprintf('%d', neuroelf_build);
np = neuroelf_path('cache');
ny = datestr(now, 10);
if exist([np filesep '.RCV'], 'file') > 0
    isrc = deblank(asciiread([np filesep '.RCV']));
else
    isrc = '';
end

% text design
lsize = 108;
ldesign = struct('bcolor', [0, 0, 0], 'color', [216, 128, 255], 'color2', [96, 96, 208], ...
    'colorgrad', [0.5, 0, 0.5, 1], 'shadow', [6, 9, 6, 0.5], 'padding', 40, ...
    'emboss', [4, 0.8], 'glow', [12, 3, 0.05], 'gcolor', [64, 32, 192], 'gcolor2', [128, 32, 96]);
fsize = 40;
design = struct('bcolor', [255, 255, 255], 'color', [64, 32, 255], 'color2', [32, 16, 96], ...
    'colorgrad', [0.5, 0, 0.5, 1], 'shadow', [4, 6, 6, 0.666], 'padding', 24, ...
    'emboss', [2, -0.3]);

netext = image_font(['NeuroElf v' nv isrc ' (build ' nb ')'], 'Harrington', round(1.6*fsize), design);
netext2 = image_font(['NeuroElf v' nv isrc ' (build ' nb ')'], 'Harrington', round(1.2*fsize), design);
ct21 = ceil(0.5 * size(netext2, 1));
ct2 = ct21 + findfirst(all(all(netext2(ct21:end, :, :) == 255, 3), 2));
netext2(ct2:end, :, :) = [];

% helpers
htext = {' '; ' '; ' '; ...
    ['(c) 2006 - ' ny ' Jochen Weber, Columbia University in the City of New York']; ' '; ...
    ['with the help of Aya Ben-Yakov, Robert Bittner, Ben Bowles, Hester Breman, ', ...
     'Shuquan Chen, Barry Cohen, Noga Cohen, Jason Craggs, Juliet']; ...
    ['Davidow, Cameron DeLeone, Federico DeMartino, Bryan Denny, Bruce Doré, ', ...
     'Fabrizio Esposito, Owen Footer, Elia Formisano, Michael Gilead,']; ...
    ['Ro''ee Gilron, Rainer Goebel, Kim Goodyear, Armin Heinecke, Chelsea Helion, ', ...
     'Katie Insel, Nir Jacoby, Zoran Josipovic, Igor Kagan, Hedy Kober,']; ...
    ['Jan Willem Koten, Ethan Kross, Frank Krueger, Ifat Levy, Rebecca Martin, ', ...
     'Maggie Mae Mell, Carmen Morawetz, Nasir Naqvi, Erik Nook,']; ...
    ['Niv Noy, Kevin Ochsner, Pim Pullens, John Pyles, Song Qi, Jenna Reinen, ', ...
     'Alard Roebroeck, Juan Sanchez, Ajay Satpute, Jen Silvers,']; ...
    ['Erez Simony, Jared van Snellenberg, Tor Wager, Noam Zerubavel, Bin Zhang, ', ...
     'and all users who report bugs and suggest new features!']; ...
    ' '; ...
    'For additional information, please visit: http://neuroelf.net/'};
fsize = repmat(fsize(1), size(htext));
fsize([2, 3]) = 0.75 * fsize(1);
fsize(5:12) = 0.6 * fsize(1);

% external support
htext2 = {' '; ' '; ...
    'I''d like to thank Xilin Shen, Emily Finn, and Monica Rosenberg for granting me permission'; ...
    'to incorporate their 268-parcel connectivity atlas into NeuroElf'; ...
    '(for more information, please see https://www.nitrc.org/frs/?group_id=51)'; ' '; ...
    'Thanks go to Luke Chang and Tor Wager for allowing their PINES pattern to be used'; ...
    '(see http://neurovault.org/collections/503/ and http://lukejchang.com/)'; ' '; ...
    'I would also like to extend my thanks to Chris Rorden for giving his blessing to bundle'; ...
    'his DICOM to NIFTI converter (dcm2nii, which is part of MRIcron) with NeuroElf!'; ...
    '(see http://people.cas.sc.edu/rorden/mricron/dcm2nii.html for details)'};
fsize2 = repmat(0.9 * fsize(1), size(htext2));
fsize2([6, 9]) = 0.4 * fsize(1);
htext3 = {' '; ' '; ...
    'A big thank-you goes to Chih-Jen Lin for allowing the re-use of his libSVM code'; ...
    '(for more information, see http://www.csie.ntu.edu.tw/~cjlin/libsvm/)'; ' '; ...
    'Also, thanks go to Aapo Hyvärinen for letting me integrate the FastICA algorithm into the'; ...
    'toolbox (please see http://research.ics.aalto.fi/ica/fastica/)'; ' '; ...
    'Finally, thanks go to Dirk-Jan Kroon, for permitting to re-use his rendering code'; ...
    '(available at http://mathworks.com/matlabcentral/fileexchange/21993-viewer3d)'; ' '; ...
    'And naturally to all users of NeuroElf who report bugs and suggest features!'};
fsize3 = repmat(0.9 * fsize(1), size(htext3));
fsize3(1:2) = 1.05 * fsize(1);
fsize3([5, 8, 11]) = 0.4 * fsize(1);

% full text
[stext1, stext1a] = image_font(['NeuroElf v' nv isrc], 'Harrington', lsize, ldesign);
stext1a = uint8(round(224 .* stext1a));
stext2 = image_font(htext, '', fsize, design);
stext3 = image_font(htext2, '', fsize2, design);
stext4 = image_font(htext3, '', fsize3, design);

% patch
stext2(1:size(netext, 1), 1:size(netext, 2), :) = netext;
stext3(1:size(netext2, 1), 1:size(netext2, 2), :) = netext2;
stext4(1:size(netext2, 1), 1:size(netext2, 2), :) = netext2;

% enforce size
if size(stext1, 1) > isize(2)
    stext1(isize(2)+1:end, :, :) = [];
    stext1a(isize(2)+1:end, :) = [];
elseif size(stext1, 1) < isize(2)
    stext1(end+1:isize(2), :, :) = 0;
    stext1a(end+1:isize(2), :) = 0;
end
if size(stext1, 2) > isize(1)
    stext1(:, isize(1)+1:end, :) = [];
    stext1a(:, isize(1)+1:end) = [];
elseif size(stext1, 2) < isize(1)
    stext1(:, end+1:isize(1), :) = 0;
    stext1a(:, end+1:isize(1)) = 0;
end
if size(stext2, 1) > isize(2)
    stext2(isize(2)+1:end, :, :) = [];
elseif size(stext2, 1) < isize(2)
    stext2(end+1:isize(2), :, :) = 255;
end
if size(stext2, 2) > isize(1)
    stext2(:, isize(1)+1:end, :) = [];
elseif size(stext2, 2) < isize(1)
    stext2(:, end+1:isize(1), :) = 255;
end
if size(stext3, 1) > isize(2)
    stext3(isize(2)+1:end, :, :) = [];
elseif size(stext3, 1) < isize(2)
    stext3(end+1:isize(2), :, :) = 255;
end
if size(stext3, 2) > isize(1)
    stext3(:, isize(1)+1:end, :) = [];
elseif size(stext3, 2) < isize(1)
    stext3(:, end+1:isize(1), :) = 255;
end
if size(stext4, 1) > isize(2)
    stext4(isize(2)+1:end, :, :) = [];
elseif size(stext4, 1) < isize(2)
    stext4(end+1:isize(2), :, :) = 255;
end
if size(stext4, 2) > isize(1)
    stext4(:, isize(1)+1:end, :) = [];
elseif size(stext4, 2) < isize(1)
    stext4(:, end+1:isize(1), :) = 255;
end

% save
sp = neuroelf_path('splash');
imwrite(stext1, [sp filesep 'text1.png'], 'Alpha', stext1a);
imwrite(stext2, [sp filesep 'text2.jpg'], 'Quality', 75);
imwrite(stext3, [sp filesep 'text3.jpg'], 'Quality', 75);
imwrite(stext4, [sp filesep 'text4.jpg'], 'Quality', 75);
