neuroelf_setup
clear classes;
clc
neuroelf_gui
neuroelf_gui('openfile', '/Volumes/mrp/Imaging/Jochen/MRP_16subs_average_ICBM.vmr');
neuroelf_gui('openfile', '/Volumes/mrp/Imaging/Jochen/MRP01/MRP01_RUN2_MNI_SEG.vtc');
neuroelf_gui('openfile', '/Volumes/mrp/Jason/Jochen/MRP_16subs_Jason_and_Jochen.vmp');
neuroelf_gui('openfile', '/Volumes/mrp/Imaging/Jochen/group/MRP_16subs_OLS+ROB_SEG_redone_cue_block_rating.vmp');
neuroelf_gui('openfile', '/Volumes/mrp/Imaging/Jochen/group/MRP_16subs_ROB_SEG_redone_cue_block_rating.glm');
neuroelf_gui('remote', 'listen');
