function normalize(input_file,output_file,thresh)
% NORMALIZE(INPUT_FILE, OUTPUT_FILE [, THRESH])
%
% OUTPUT_FILE is INPUT_FILE divided by the average of INPUT_FILE above
% THRESH. If THRESH is not supplied, uses half the max of INPUT_FILE. 
d=fmris_read_image(deblank(input_file));
if nargin<=2
   thresh=max(d.data(:))/2
end
av=mean(d.data(find(d.data>thresh)))
d.data=d.data/av;
d.file_name=deblank(output_file);
fmris_write_image(d);
return

