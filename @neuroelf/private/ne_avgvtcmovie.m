function varargout = ne_avgvtcmovie(varargin)

% structure:
% - inquire about either voxel-space (VMR), surface-space (SRF), or render
% - also show time course plot (with vertical indicator) in movie?
% - join VTCs if needed (into single AvgVTC)
% - select conditions to show
% - sampling window (and advancing speed as seconds/s)
% - thresholding options
% - sample the VTC(s) if needed -> MTC
% - if voxel space: three-slice (orthogonal) projection or multi-slice view
% - if surface space: which surfaces and angles
% - if multiple conditions (time courses): stacking or co-coloring
% - video options (frame rate, default: 24/s)
% - size options (1080full, 800x600, etc.)
% - annotation (incl. time axis) in frames
% - target filename (of movie file)
% - advanced option: allow to select image to overlay/underlay within
%   movie frame for specific times (e.g. experimental slides, like images
%   rating scale, fixation cross, etc.)
% - advanced option: allow to create flexible time axis (with breaks,
%   reversals, etc.)

% processing, frame-by-frame screenshoting, imread, addframe, 
