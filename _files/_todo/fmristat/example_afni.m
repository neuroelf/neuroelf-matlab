% Looking at the fMRI data using pca_image

input_file='c:/keith/kwafnidata/r1_time+orig.BRIK';
mask_file='c:/keith/kwafnidata/r1_time+orig.BRIK';
mask_thresh=fmri_mask_thresh(mask_file);
pca_image(input_file, [], 4, mask_file, mask_thresh);
saveas(gcf,'c:/keith/test/figs_afni/figpca1.jpg');

exclude=[1 2 3 4];
pca_image(input_file, exclude, 4, mask_file);
saveas(gcf,'c:/keith/test/figs_afni/figpca2.jpg');

X_remove=[ones(68,1) (1:68)'];
X_interest=[zeros(4,8); repmat(eye(8),8,1)];
pca_image(input_file, exclude, 4, mask_file, [], [], [], X_remove, X_interest);
saveas(gcf,'c:/keith/test/figs_afni/figpca3.jpg');

% Making the design matrices using fmridesign

[err, Info]=BrikInfo('c:/keith/kwafnidata/r1_time+orig.HEAD');
Info
TR=5;
frametimes=(0:67)*TR;
slicetimes=Info.TAXIS_OFFSETS/1000;
eventid=repmat(1,8,1);
eventimes=((0:7)'*8+4)*TR;
duration=ones(8,1)*4*TR;
height=ones(8,1);
events=[eventid eventimes duration height] 
X_cache=fmridesign(frametimes,slicetimes,events);
clf;
plot(squeeze(X_cache.X(:,:,1,4)),'LineWidth',2)
xlabel('frame number')
ylabel('response')
saveas(gcf,'c:/keith/test/figs_afni/figdesign.jpg');

% Analysing one run with fmrilm

contrast=[1];
exclude=[1 2 3 4];
which_stats='_mag_t _mag_ef _mag_sd _cor _fwhm';

input_file='c:/keith/kwafnidata/r1_time+orig.BRIK';
output_file_base='c:/keith/test/results_afni/r1_time';

fmrilm(input_file, output_file_base, X_cache, contrast, exclude, which_stats);

% Visualizing the results using view_slices, glass_brain and blob_brain

t_file='c:/keith/test/results_afni/r1_time_mag_t+orig.BRIK';
m=fmris_read_image(t_file,8,1);
clf;
imagesc(m.data',[-6 6]); colorbar; axis xy; colormap(spectral);
saveas(gcf,'c:/keith/test/figs_afni/figslice.jpg');

mask_file=input_file;
clf;
view_slices(t_file,mask_file,[],7,1,[-6 6]);
saveas(gcf,'c:/keith/test/figs_afni/figviewslice.jpg');

clf;
view_slices(t_file,mask_file,[],0:15,1,[-6 6]);
saveas(gcf,'c:/keith/test/figs_afni/figviewslices.jpg');

glass_brain(t_file,3,mask_file);
saveas(gcf,'c:/keith/test/figs_afni/figlassbrain.jpg');

clf;
blob_brain(t_file,5,mask_file);
title('T>5, coloured by effect (%BOLD)');
saveas(gcf,'c:/keith/test/figs_afni/figblobrain.jpg');

cor_file='c:/keith/test/results_afni/r1_time_cor+orig.BRIK';
clf;
view_slices(cor_file,mask_file,0,7,1,[-0.15 0.35]); 
saveas(gcf,'c:/keith/test/figs_afni/figcor.jpg');

fwhm_file='c:/keith/test/results_afni/r1_time_fwhm+orig.BRIK';
clf
view_slices(fwhm_file,mask_file,0,7,1,[0 20]);
saveas(gcf,'c:/keith/test/figs_afni/figfwhm.jpg');

% F-tests

contrast.T=[0 1 0 0;
            0 0 1 0;
            0 0 0 1]
which_stats='_mag_F';
output_file_base='c:/keith/test/results_afni/r1_time_drift';
fmrilm(input_file,output_file_base,X_cache,contrast,exclude,which_stats,cor_file);
clf;
view_slices('c:/keith/test/results_afni/r1_time_drift_mag_F+orig.BRIK',mask_file,[],0:15,1,[0 50]);
saveas(gcf,'c:/keith/test/figs_afni/figdrift.jpg');

% Thresholding the tstat image with stat_threshold and fdr_threshold

[search_volume, num_voxels]=mask_vol(mask_file)
stat_threshold(search_volume,num_voxels,6.6393,52)
fdr_threshold(t_file,[],mask_file,[],52)

% Producing an SPM style summary with stat_summary

stat_summary(t_file, fwhm_file, [], mask_file);
saveas(gcf,'c:/keith/test/figs_afni/figlassbrainmulti.jpg');

% Confidence regions for the spatial location of local maxima using conf_region

conf_region(t_file,4.97,mask_file)
conf_file='c:/keith/test/results_afni/r1_time_mag_t_95cr+orig.BRIK';
clf;
view_slices(conf_file,mask_file,[],0:15)
clf;
blob_brain(conf_file,4.97,conf_file,4.97)
title('Approx 95% confidence regions for peak location, coloured by peak height')
saveas(gcf,'c:/keith/test/figs_afni/figconfregion.jpg');

% Extracting values from a minc file using extract

voxel=[32  26   7]
ef=extract(voxel,'c:/keith/test/results_afni/r1_time_mag_ef+orig.BRIK')  
sd=extract(voxel,'c:/keith/test/results_afni/r1_time_mag_sd+orig.BRIK')
ef/sd

[df, spatial_av]=fmrilm(input_file,[],[],[],exclude);
ref_data=squeeze(extract(voxel,input_file))./spatial_av*100;
fitted=mean(ref_data)+ef*X_cache.X(:,1,1,voxel(3)+1);
clf;
plot(frametimes,[ref_data fitted],'LineWidth',2); 
xlim([15 max(frametimes)]);
legend('Reference data','Fitted values');
xlabel('time (seconds)');
ylabel('fMRI response, percent');
title(['Observed (reference) and fitted data, ignoring trends, at voxel ' num2str(voxel)]);
saveas(gcf,'c:/keith/test/figs_afni/figfit.jpg');

% Estimating the time course of the response

eventid=kron(ones(8,1),(1:8)');
eventimes=((1:64)'+3)*TR;
duration=ones(64,1)*TR;
height=ones(64,1);
events=[eventid eventimes duration height]
X_bases=fmridesign(frametimes,slicetimes,events,[],zeros(1,5));
contrast=[eye(8)-ones(8)/8];
num2str(round(contrast*100)/100)
which_stats='_mag_ef _mag_sd _mag_F';
output_file_base=['c:/keith/test/results_afni/r1_time01';
'c:/keith/test/results_afni/r1_time02';
'c:/keith/test/results_afni/r1_time03';
'c:/keith/test/results_afni/r1_time04';
'c:/keith/test/results_afni/r1_time05';
'c:/keith/test/results_afni/r1_time06';
'c:/keith/test/results_afni/r1_time07';
'c:/keith/test/results_afni/r1_time08'];
df=fmrilm(input_file,output_file_base,X_bases,contrast,exclude,which_stats,cor_file)

stat_threshold(search_volume,num_voxels,6,df.F)
lm=locmax('c:/keith/test/results_afni/r1_time01_mag_F+orig.BRIK',7.38);
num2str(lm)

values=extract(voxel,output_file_base,'_mag_ef+orig.BRIK')
sd=extract(voxel,output_file_base,'_mag_sd+orig.BRIK')

time=(1:(8*TR*10))/10;
X_hrf=fmridesign(time,0,[1 0 20 1]);
hrf=squeeze(X_hrf.X(:,1,1));
plot((0:8)*TR+slicetimes(voxel(3)+1),values([1:8 1]),'k', ...
   [0:7; 0:7]*TR+slicetimes(voxel(3)+1), [values+sd; values-sd],'g', ...
   time,[ones(1,200) zeros(1,200)],'r', ...
   time,hrf*ef,'g','LineWidth',2);
legend('Estimated response','Modeled response');
xlabel('time (seconds) from start of epoch');
ylabel('fMRI response, percent');
title(['Estimated and modeled response at voxel ' num2str(voxel)]);
saveas(gcf,'c:/keith/test/figs_afni/figmodelresp.jpg');

% Estimating the delay 

contrast=[1]
which_stats='_del_t _del_ef _del_sd'
output_file_base=['c:/keith/test/results_afni/r1_time']
fmrilm(input_file, output_file_base, X_cache, contrast, exclude, which_stats, cor_file);

slice=7; 
subplot(2,2,1);
view_slices('c:/keith/test/results_afni/r1_time_mag_t+orig.BRIK',mask_file,[],slice,1,[-6 6]);
subplot(2,2,2);
view_slices('c:/keith/test/results_afni/r1_time_del_ef+orig.BRIK',mask_file,[],slice,1,[-3 3]);
subplot(2,2,3);
view_slices('c:/keith/test/results_afni/r1_time_del_sd+orig.BRIK',mask_file,[],slice,1,[0 6]);
subplot(2,2,4);
view_slices('c:/keith/test/results_afni/r1_time_del_t+orig.BRIK',mask_file,[],slice,1,[-6 6]);
saveas(gcf,'c:/keith/test/figs_afni/figdelay.jpg');

clf;
blob_brain('c:/keith/test/results_afni/r1_time_mag_t+orig.BRIK',5, ...
   'c:/keith/test/results_afni/r1_time_del_ef+orig.BRIK');
title('Delay (secs) of hot stimulus where T > 5')
saveas(gcf,'c:/keith/test/figs_afni/figdelay3D.jpg');

% Effective connectivity of all voxels with a reference voxel 

output_file_base='c:/keith/test/results_afni/r1_time_connect';
contrast.C=1;
which_stats='_mag_t _mag_ef _mag_sd'
ref_times=frametimes'+slicetimes(voxel(3)+1);
confounds=fmri_interp(ref_times,ref_data,frametimes,slicetimes);
fmrilm(input_file, output_file_base, X_cache, contrast, exclude, which_stats, cor_file, [], confounds);
clf;
view_slices('c:/keith/test/results_afni/r1_time_connect_mag_t+orig.BRIK',mask_file,[],0:15,1,[-6 6]);
saveas(gcf,'c:/keith/test/figs_afni/figconnect.jpg');

output_file_base='c:/keith/test/results_afni/r1_time_connect_stim';
contrast.C=[0 1];
X_confounds=XC(X_cache,confounds);
fmrilm(input_file, output_file_base, X_cache, contrast, exclude, ...
    which_stats, cor_file, [], X_confounds);
clf;
view_slices('c:/keith/test/results_afni/r1_time_connect_stim_mag_t+orig.BRIK',mask_file,[],0:15,1,[-6 6]);
saveas(gcf,'c:/keith/test/figs_afni/figconnect_stim.jpg');

% Higher order autoregressive models 

which_stats='_mag_t _mag_ef _mag_sd _cor _AR'
contrast=[1];
output_file_base='c:/keith/test/results_afni/r1_time_ar4';
fmrilm(input_file, output_file_base, X_cache, contrast, exclude, ...
    which_stats, [], [], [], [], [], [], 4);
clf;
view_slices('c:/keith/test/results_afni/r1_time_ar4_AR+orig.BRIK',mask_file,0,7,1:4,[-0.15 0.35]); 
saveas(gcf,'c:/keith/test/figs_afni/figar4.jpg');

