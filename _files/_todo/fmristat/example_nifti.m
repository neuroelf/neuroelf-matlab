input_file='c:/keith/nifti/nifti/filtered_func_data.nii'
fmris_read_image(input_file,0,0)

view_slices(input_file,input_file,1000,0:20,1)
pca_image(input_file,[],4,input_file,1000)

frametimes=(0:179)*3;
slicetimes=(0:20)/21*3;
events=[ones(9,1) (1:9)'*20-10 ones(9,1)*10 ones(9,1)]
X_cache=fmridesign(frametimes,slicetimes,events);

contrast=[1];
exclude=[];
which_stats='_mag_t _mag_sd _mag_ef _cor _fwhm';
output_file_base=['c:/keith/nifti/nifti/test']
fmrilm(input_file, output_file_base, X_cache, contrast, exclude, which_stats);

fmris_read_image([output_file_base '_mag_t.nii'],0,0)
clf; view_slices([output_file_base '_mag_t.nii'])
clf; view_slices([output_file_base '_mag_sd.nii'])
clf; view_slices([output_file_base '_mag_ef.nii'])
clf; view_slices([output_file_base '_cor.nii'])
clf; view_slices([output_file_base '_fwhm.nii'],[],[],0:20,1,[0 20])

d=fmris_read_image('c:/keith/nifti/nifti/zstat.nii')
d.file_name='c:/keith/nifti/nifti/zstat2.nii'
d.fwhm=[11 35]
fmris_write_image(d)
fmris_read_image('c:/keith/nifti/nifti/zstat2.nii',0,0)
view_slices('c:/keith/nifti/nifti/zstat2.nii')

for i=1:21
d=fmris_read_image('c:/keith/nifti/nifti/zstat.nii',i,1)
d.file_name='c:/keith/nifti/nifti/zstat2.nii'
d.fwhm=[11 35]
fmris_write_image(d,i,1)
end
fmris_read_image('c:/keith/nifti/nifti/zstat2.nii',0,0)
view_slices('c:/keith/nifti/nifti/zstat2.nii')

file='c:/keith/nifti/nifti/filtered_func_data.nii'



