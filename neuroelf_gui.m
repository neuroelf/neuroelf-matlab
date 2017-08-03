%NEUROELF_GUI  NeuroElf graphical user interface (GUI).
%
%   To invoke the GUI, simply enter neuroelf_gui without any arguments.
%
%   In addition, you can call (and thus script) almost any action the GUI
%   can perform by adding an 'action' argument (followed by any additional
%   optional arguments) to the call:
%
%   FORMAT: [out, ...] = neuroelf_gui([action, [ arguments]])
%
%   To find out more about an action, please try
%
%       action      string representing an action to perform
%       arguments   additional arguments provided to action function
%
% Output fields:
%
%       out, ...    any output provided by action function
%
% Note: if called with arguments, neuroelf_gui tries to resolve a function
%       with the location NEUROELF_FOLDER/@neuroelf/private/ne_ACTION.m
%
% For help on available sub-functions, please use
%
% neuroelf_gui('help') and neuroelf_gui('help', 'action')
%
% Matlab toolbox developed, written and maintained by Jochen Weber
% with the help and inspiration of (affiliations may be outdated!)
%
% - Aya Ben-Yakov, Weizmann Institute of Science, Rehovot, Isreal
% - Robert Bittner, Goethe University, Frankfurt aM, Germany
% - Ben Bowles, UC Berkeley, CA, USA
% - Hester Breman, Brain Innovation, B.V., Maastricht, Netherlands
% - Barry Cohen, NYU, New York, NY, USA
% - Noga Cohen, Columbia University, New York, NY, USA
% - Jason Craggs, University of Missouri, Columbia, MO, USA
% - Juliet Davidow, Harvard University, Cambridge, MA, USA
% - Cameron DeLeone, Yale University, New Haven, CT, USA
% - Federico DeMartino, University of Maastricht, Netherlands
% - Bryan Denny, Rice University, Houston, TX, USA
% - Bruce Doré, Columbia University, New York, NY, USA
% - Fabrizio Esposito, University of Salerno, Fisciano, SA, Italy
% - Owen Footer, Columbia University, New York, NY, USA
% - Elia Formisano, University of Maastricht, Netherlands
% - Michael Gilead, Ben-Gurion University, Beersheba, Israel
% - Ro'ee Gilron, Tel Aviv University, Tel Aviv, Israel
% - Rainer Goebel, Brain Innovation, B.V., Maastricht, Netherlands
% - Kim Goodyear, George Mason University, Fairfax, VA, USA
% - Armin Heinecke, Brain Innovation, B.V., Maastricht, Netherlands
% - Chelsea Helion, Columbia University, New York, NY, USA
% - Katie Insel, Harvard University, Cambridge, MA, USA
% - Nir Jacoby, Columbia University, New York, NY, USA
% - Zoran Josipovic, NYU, New York, NY, USA
% - Igor Kagan, Deutsches Primatenzentrum, Goettingen, Germany
% - Hedy Kober, Yale University, New Haven, CT, USA
% - Jan Willem Koten, RWTH Aachen, Germany
% - Ethan Kross, University of Michigan, Ann Arbor, MI, USA
% - Frank Krueger, George Mason University, Fairfax, VA, USA
% - Ifat Levy, Yale University, New Haven, CT, USA
% - Rebecca Martin, Columbia University, New York, NY, USA
% - Maggie Mae Mell, Yale University, New Haven, CT, USA
% - Carmen Morawetz, FU Berlin, Germany
% - Nasir Naqvi, NY State Psychiatric Institute, New York, NY, USA
% - Erik Nook, Harvard University, Cambridge, MA, USA
% - Niv Noy, Weizmann Institute of Science, Rehovot, Isreal
% - Kevin Ochsner, Columbia University, New York, NY, USA
% - Pim Pullens, Antwerp University, Netherlands
% - John Pyles, Carnegie Mellon University, Pittsburgh, PA, USA
% - Song Qi, Columbia University, New York, NY, USA
% - Jenna Reinen, Yale University, New Haven, CT, USA
% - Alard Roebroeck, University of Maastricht, Netherlands
% - Juan Sanchez, NY State Psychiatric Institute, New York, NY, USA
% - Ajay Satpute, Pomona College, Claremont, CA, USA
% - Jen Silvers, UCLA, Los Angeles, CA, USA
% - Erez Simony, Princeton University, Princeton, NJ, USA
% - Jared van Snellenberg, Columbia University, New York, NY, USA
% - Tor Wager, University of Colorado, Boulder, CO, USA
% - Noam Zerubavel, Columbia University, New York, NY, USA
% - Bin Zhang, NYU, New York, NY, USA
%
% additional thanks go to the SCAN Unit of Columbia University to allow
% me to continue working on this toolbox, and to all users who report
% bugs and request features.
%
%
% Menu commands:
%
% - File
%   - Open               load any supported document, anatomical and
%                        statistical projects will be added to one of two
%                        lists of available objects, from which to choose
%   - Open file as stats allows to load general files (HEAD/HDR/NII) as
%                        a statistics (instead of as an anatomical) file
%
%   - Recently loaded    list with recently loaded objects (in four
%     ... objects        different categories: slicing, stats, surface,
%                        and surface stats)
%
%   - Colin-27 dataset   various objects that are created by calling
%                        neuroelf_makefiles (possibly with an 'all' input)
%
%   - NeuroSynth maps    maps downloadable (and previously downloaded) from
%                        Tal Yarkoni's http://neurosynth.org/ (please note
%                        that for now, this downloads OLDER maps, not the
%                        most up-to-date version!!)
%
%   - Clone slicing      create a copy of the current anatomical dataset
%     object             e.g. for attempting to draw (the original dataset
%                        could then be set as an "underlay" object for
%                        additional guidance)
%   - New slicing object creates a new (empty) VMR, e.g. for drawing
%   - Reload slicing     reloads the currently selected anatomical dataset
%     object             from disk (e.g. when too much smoothing was done)
%   - Set underlay       allows to select another anatomical dataset to be
%     object             displayed below (mix into) the slicing
%
%   - Clone stats object same as with the anatomical (testing operations)
%   - Reload stats object equally, reloading the stats object (undo)
%
%   - Save/as            saves the currently selected anatomical dataset
%                        (with Save as... asking for a new filename)
%
%   - Save stats/as      saves the currently selected statistics dataset
%                        (with Save as... asking for a new filename)
%
%   - Save text output   allows to save the text shown in the lower left
%     (tables/betas)     edit field to a text file (without copy/paste)
%
%   - Reload scenery     reload a previously created scenery (MAT) file
%     file               containing surfaces, surface statistics, and the
%                        corresponding configurations
%
%   - Import Analyze to  import a structural dataset from a HDR/NII (or
%     VMR                AFNI HEAD/BRIK) file into VMR format
%   - Import SPM maps to import several SPM-based maps (HDR/NII files)
%     VMP                into a single VMP (for convenient storage, etc.)
%   - Import RFX-GLM     ability to convert an SPM-based analysis stream
%     from SPM.mat files of 1st-level regressions into a NeuroElf GLM file
%
%   - Close all files    cloes all open objects (i.e. start new session)
%
%   - Options
%     - ClusterTable     settings related to cluster table generation
%     - Draw on event    when to actually draw (e.g. mouse, keyboard, etc.)
%     - NeuroElf listener remote configuration (requires HTTP server!)
%     - Orientation      switch between neurological and radiological
%     - PLP-Lookup       settings for how to lookup from current PLP
%     - Renderer         choose bewteen OpenGL and zbuffer renderer
%     - Spatial          allows to attach an SPM-based *_sn.mat file to
%       normalization    either anatomical or statistical dataset for
%                        on-the-fly spatial normalization (slice display)
%     - Surface view     set options related to Surface view (and also
%       options          for reconstruction of surfaces from VMRs)
%
%     - Grayscale colors set a different LUT (color table) for anatomical
%                        datasets (e.g. for border/inhomogeneity detection)
%     - Stats colors     load a different LUT (for the display of non-VMP
%                        statistical maps, as well as to set the .LUTName
%                        field of the currently selected VMP.Map)
%     - Edit LUT         manually choose colors (RGB editor: colorpicker)
%
%     - Underlay         choose how the main anatomical and underlay
%       blending mode    anatomical datasets are blended
%
%     - Compute inst.    turn on/off the computation of seed correlations
%       seed correlation at the current position (for VTCs only!)
%     - Create new VMR in  option as to whether File -> New VMR creates the
%       0.5mm resolution   dataset in 1mm or 0.5mm resolution
%     - Extended map     if checked, VMP (and other) map names will be
%       names in list    extended by type, d.f. information in the list
%     - Linked browsing  if checked, setting the cursor position will be
%       with satellites  synchronized across all NeuroElf UI windows
%     - Maximum distant  when displaying exactly two (2) statistical maps
%       color for two    the colors for voxels that surpass the threshold
%       joined maps      in both maps will be re-computed in HSV space
%     - Show stats       display a (set of) heat-color bar(s) next to a
%       thresh bars      dataset on which statistic maps are overlaid
%     - Echo function    print to the console important calls which can be
%       calls/methods    used as templates for user-defined scripts
%
%   - Exit               close the GUI
%
% - Object-related menus between the File and Analysis menus, NeuroElf
%                        displays a menu relevant for the currently
%                        selected dataset (e.g. VMR, VMP, SRF, SMP, etc.)
%
% - VMR                  optionally visible menu for VMR datasets
%   - Border reconstruction creates a SRF object representing the border
%   - Correct inhomogeneity corrects the inhomogeneity (requires V16 data!)
%   - Recompute 8bit from 16bit takes existing V16 data to recompute VMR
%
%   - Export to Nifti    create a NII file from current VMR
%   - Export as RGB Nifti creates a NII file with 3 values per voxel (RGB)
%
% - GLM
%
% - Analysis
%   - Beta plot          for GLMs, bring up additional window that plots
%                        beta weights of a GLM at the cursor position
%   - Contrast manager   compute t-statistic maps
%
% - Visualization
%   - Create montage     create series of slices through the currently
%                        selected slicing object (incl. stats maps/options)
%   - Create surface     take a snapshot from the surface window (axes),
%     snapshot           also incl. the stats maps/options
%
% - Tools
%   - SPM -> BV conversion
%     - Import VMR       converts one anatomical image (img/nii) into a VMR
%     - Import stats     converts one or several spmT/spmF images into a
%     VMP
%     - Import SPM.mats  converts a list of 1st-level SPM.mat files (after
%                        estimation) into a random-effects GLM file
%     - SPM.mat -> PRT   extract onsets from an SPM.mat file and create one
%                        or several PRT(s)
%     - SPM.mat -> SDM   extract the design matrices of the runs and create
%                        the corresponding SDM files
%
%   - alphasim           run a MC-simulation on random maps (smoothed) to
%                        determine the cluster threshold for several
%                        uncorrected thresholds
%   - fmriquality        assess the quality of an fMRI run
%   - renamedicom        flexibly rename several DICOM files (supports
%                        subdirectories)
%   - SPM preprocessing  configure (and run) the preprocessing on several
%                        subjects' datasets
%   - tdclient           manually lookup coordinates to labels, cubes, etc.
%
% - Help
%   - About              bring up this dialog
%
%
%
% Interactive commands (clicks and edits):
%
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % (1) anatomical %^             %             %
%   % (2) statistics %7      (8)    %     (9)    v%(15)
%   %                %v   sagittal  %   coronal  i%
%   % (3) stats maps %     slicing  %   slicing  e%
%   %                %d             %            w%
%   %                %r  %%%%%%%%%%%%%%%%%%%%%%%  %
%   % (4) clusters   %a     (11)    %            c%(16)
%   %                %w  coordinates%    (10)    o%
%   %                %i  %%%%%%%%%%%%    axial   n%
%   % (5) tables and %n     (12)    %   slicing  t%
%   %     extracted  %g thresholding%             %
%   %     betas      %  (13) sampled values, info %
%   %                %                            %
%   % (6) progress   %  (14) time course display  %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% - selecting the anatomical object (1), which includes 4D HDRs and VTCs,
%   shows the slicing of that object; for functional datafiles (VTCs, HDRs,
%   and FMRs), the time course at the cursor is also displayed (14)
% - the orientation button allows to specify object-specific settings
% - clicking the "x" button will remove this object from the UI (and, if
%   it was loaded by the GUI, the contents will be cleared from memory)
%
% - selecting the statistical object (2) shows a list of maps (3) available
%   in that object for display, and makes certain buttons and options
%   selectable, to the right of to (2) and (3)
% - the orientation button works in parallel to the anatomical object
% - and so does the "x" button for closing objects
% - if a statistical object contains references to VTCs (for now only GLMs)
%   those can be selected and will then be loaded (as anatomical datasets)
% - the additional buttons can be used to change whether the expected
%   response will also be plotted in the time course portion (14)
% - and whether or not to remove nuisance variance from the time courses
%
% - selecting a single available map (3) enables the threshold controls (12)
% - if multiple maps are selected (COMMAND + click), several maps
%   are being overlaid (with either mixing the colors or such that maps
%   occuring later in the list are having precendence over previous maps);
%   for multiple selected maps, thresholding controls (12) are disabled
% - the up and down buttons re-order the maps in the VMP/HEAD file
% - the properties button allows to change some settings of the (single)
%   selected map (in a VMP/SMP/HEAD container)
% - the trashcan button removes the selected maps from the VMP/SMP/HEAD
% - the function ("f") button allows to compute a mathematical expression
%   over the selection ($INDEX) and all available (#INDEX) maps, e.g.
%   $1 .* (abs($2) > 1.96) would create a new map from the first selected
%   map masked to those voxels that reach abs(stat)>1.96 in the second map
% - the surface-brain button samples the currently selected maps into
%   surface space (with the currently loaded SRF)
% - the beta-plotter button opens the GLM beta plotter; if the VMP map
%   contains information about the contrast, those settings will be applied
% - the ellipse-button propagates the selection of a map across subjects
%   of a RFX-GLM
%
% - selecting a single cluster (4) will set the current position to the
%   first (peak) voxel of a (sub-) cluster; if enabled, it will also
%   extract the betas from the currently selected GLM
% - using the marking button, voxels of the selected clusters will be
%   drawn into the VMR, with configurable color code and extent setting
% - using the sphere button, clusters can be restricted to a spheroid; also
%   if no cluster is selected, a spherical ROI can be created at the cursor
% - clicking the ATL (atlas) button allows to add ROIs from the Talairach
%   atlas, converted to ICBM if the checkbox in (11) is selected
% - using the crosshair button will set the cursor to the closest peak in
%   the table (one of the first voxels across all ROIs)
% - the lens buttons zooms the slicing view to the peak of a cluster
% - the trashcan button removes the selected clusters from the VOI
% - the folder and disk buttons allow loading and saving of the VOI/clusters
% - clicking the table-button next to the cluster list extracts betas
%   from the selected clusters (a GLM must be selected)
%
% - the edit field (5) can be used to select and copy text (i.e. the output
%   of cluster table and beta extraction operations)
%
% - the progress bar (6) will be made visible during long-running tasks
%   (e.g. contrast computations, alphasim, etc.), and in addition the
%   UI title bar will also indicate the progress in such situations
%
% - the drawing tools (7) allow to "paint" into the selected VMR/NII
% - clicking the cross-hair button disables the drawing mode (browsing)
% - the pen (or 3D pen) enables 2D (3D) drawing (configurable settings);
%   when drawing with a spherical pen, additionally the smoothness of the
%   drawing tool can be set for non-special values (e.g. to smoothen edges)
% - clicking on the bucket allows to fill a region of space starting with
%   the current position configurable over a value range
% - clicking on the border extension button will extend the painting color
%   into 1-degree neighbor voxels within a range (to grow the marking)
% - the gray-scale smoothing button can be used to smooth an anatomical
%   dataset *without* marked voxels
% - the marker-smoothing button will smooth the binary selection and then
%   reapply this to the dataset (e.g. good for smoothing edges of ROIs)
% - the within-range smoothing only smoothes values within a value range
% - the undo-pen button toggles the drawing between marker and restore mode
% - the check-mark button accepts all changes (i.e. copies the data to the
%   undo-buffer)
% - the revert button reverts all non-accepted changes (copy buffer back)
% - the reload (double revert) button reloads the content from disk
% - the green masking button reloads only voxels marked in the currently
%   configured drawing color
% - the inverse-masking button reloads only the non-marked voxels
% - the ROI button creates a new ROI in (4) from the selected voxels
% - the underlay button allows to select an underlay object (same as
%   File -> Set underlay object)
% - the V16 button toggles between displaying the regular VMRData or the
%   V16data from a VMR dataset with auto-loaded V16
%
% - clicks into either of the three panels (8, 9, or 10) sets the cursor to
%   a new position (browsing), and possibly draws at the new location
%
% - the position can also be manually set with the voxel edit boxes (11)
% - these controls also allow to set the volume index for a 4D dataset;
%   however, this is more convenient by clicking into the time course (14)
%
% - thresholds can either be set with the edit boxes or the dropdown (12),
%   which takes the number of tails into account
% - the positive/negative checkboxes toggle the corresponding stats tail
% - when the k-thresh box is checked, the displayed map will only show
%   clusters that surpass the number of contiguous voxels in size; please
%   note that this number is in voxel units, so the real-world (mm^3) size
%   of any clusters is dependent on the data resolution!
% - the alphasim button will (configurably) run the alphasim function for
%   the currently selected map; the results will be stored in the map's
%   RunTimeVars struct for later use
% - the Clustertable button will invoke the clustering of a map which then
%   produces a text table
% - if larger clusters are to be split up (by a watershed algorithm) into
%   their sub-clusters, please check the "split" box
% - whenever the icbm2tal checkbox is ticked, the UI will perform the
%   necessary coordinate conversion (using the icbm2tal and tal2icbm
%   functions of neuroelf); this applies to labels in the info box (13) as
%   well as to the extracted text tables (5) and atlas ROIs (ATL button)
% - in case you wish to retrieve the Talairach Daemon Database label for
%   peak voxels in cluster tables (4) and (5), as well as the current
%   position (11), please enable the TDclient box
% - the interpolation checkbox toggles interpolation of stats maps
% - the LUT and RGB radio buttons allow to switch between LUT and RGB
%   coloring schemes for VMP maps
% - clicking on one of the RGB color buttons will bring up a colorpicker
%   which allows to set the colors through the GUI
%
% - the view buttons (15) can be used to alter the general look and feel
%   of the GUI
% - the slicing view can be altered with selecting the combined or a
%   singular view button (right of slicing view)
%
% - the translation/rotation/scaling of the dataset (for display) can be
%   altered with the view properties button
%
% - the surface view is available via the scenery button
%
% - a click into a displayed timecourse also sets the volume number and
%   updates the display (which allows interactive volume browsing)
%
% - in the scenery listbox, all selected SRFs will be shown concurrently
%
% - if only one SRF is selected, the display properties can be set manually
%
%   the following keyboard commands are handled by the main GUI
%
%   *without* any modifier keys (neither of SHIFT, CONTROL, or ALT pressed)
%
%   'a'         - toggle current position lines on/off for slicing view
%   'b'         - switch to browse (non-drawing) mode
%   'c' / 'k'   - toggle cluster-size thresholding on/off
%   'd'         - switch to 2D drawing mode
%   'f'         - set currently selected surface to 'faces' display
%   'g'         - toggle gradient display mode on/off (for slicing view)
%   'i'         - toggle interpolation on/off (for slicing view)
%   'j'         - toggle stats-map color joining mode on/off
%   'l'         - toggle local-max splitting on/off
%   'm'         - cycle through interpolation modes ('linear', 'cubic')
%   'n'         - toggle stats alpha-thresholding setting on/off
%   'o'         - toggle orientation neurological/radiological (slicing)
%   'r'         - reset most aspects of currently displayed view
%   's'         - toggle between small/large UI size
%   't'         - for slicing view, cycle through triple- and single-slicing
%                 for surface view, toggle transparency on/off (single SRF)
%   'u'         - toggle undo-drawing mode on/off
%   'w'         - set currently selected surface to 'wireframe' display
%   'x'         - cycle cluster-extraction mode ('manual', 'single', 'multi')
%   'z'         - toggle zoom-mode on/off (slicing view)
%   '1'         - toggle positive stats tail on/off
%   '2'         - toggle negative stats tail on/off
%   '3'         - switch to 3D drawing mode
%   cursor l/r  - for slicing view, move along X-axis (left/right)
%                 for surface/render views, rotation around Z axis
%   cursor u/d  - for slicing view, move along Z-axis (up/down)
%               - for surface/render views, set zenith viewing angle
%
%   *with* SHIFT pressed
%
%   'b'         - toggle ShowThreshBars on/off
%   'i'         - toggle instantaneous seed-correlation on/off
%   'j'         - toggle max-color-dist for two stats maps on/off
%   'm'         - cycle through stats maps in current container
%   's'         - create screenshot image file
%   '1' .. '9'  - select current dataset (slicing / stats object)
%   cursor l/r  - for slicing view, move along time-axis (earlier/later)
%                 for surface/render views, translate/shift (left/right)
%   cursor u/d  - for slicing view, move along Y-axis (front/back)
%                 for surface/render views, translate/shift (up/down)
%
%   *with* CONTROL pressed
%
%   'l'         - toggle linked UIs on/off
%   'r'         - reset slicing object specific transformation matrix
%   's'         - cycle through all slicing objects
%   cursor u/d  - for surface view, increase/decrease zoom factor
%
%   *with* ALT pressed
%
%   's'         - maximize size of main GUI (factor up)
%   '1'         - switch to 3-slice (slicing) view
%   '2'         - switch to SAG single-slicing view
%   '3'         - switch to COR single-slicing view
%   '4'         - switch to TRA single-slicing view
%   '5'         - switch to surface view
%   '6'         - switch to render view (and open render UI if needed)
%   cursor u/d  - for single stats map, up or down transparency by 0.1
%
%   depending on the Operating System, the following keys will
%
%   *with* COMMAND pressed (or CONTROL for Windows)
%
%   'c'         - open the contrast manager (requires loaded GLM)
%   'd'         - open the single-level (RFX) mediation dialog
%   'e'         - open the ECG/heart/physio-data analysis dialog
%   'g'         - open the MDM::ComputeGLM dialog
%   'i'         - open the image montage creation dialog
%   'o'         - bring up the File->Open dialog
%   'p'         - open the SPM5/8-based scripted preprocessing dialog
%   'r'         - open and switch to the rendering dialog
%   's'         - save currently selected slicing var (save over!)
%   't'         - open an interactive Talairach Daemon Database client (UI)
%   'x'         - close NeuroElf GUI
%   'y'         - re-load previously stored scenery objects file
%
%
%   and the following commands are handled differently by satellite windows
%
%   'b'         - choose background color for surface view

% Using: ne_initwindow, ne_showpage.

% Version:  v1.1
% Build:    16061718
% Date:     Jun-17 2016, 6:27 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2016, Jochen Weber
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
function [varargout] = neuroelf_gui(varargin)

% global variable (used by all functions and also some xff methods)
global ne_gcfg ne_methods;

% pre-set output
varargout = cell(1, nargout);

% only if UI is not already up and running
if numel(ne_gcfg) ~= 1 || ~isstruct(ne_gcfg) || ~isfield(ne_gcfg, 'h') || ...
   ~isfield(ne_gcfg.h, 'MainFig') || ~isxfigure(ne_gcfg.h.MainFig, true)

    % deal with errors gracefully
    try

        % delete prior objects with same Tag and (re-) create and main fig
        shh = get(0, 'ShowHiddenHandles');
        set(0, 'ShowHiddenHandles', 'on');
        nemf = findobj('Tag', 'NeuroElf_MainFig');
        if ~isempty(nemf)
            set(nemf, 'CloseRequestFcn', '');
            delete(nemf);
        end
        set(0, 'ShowHiddenHandles', shh);

        % reset ne_gcfg
        ne_gcfg(:) = [];

        % show about dialog
        try
            fprintf('NeuroElf (%d): ', neuroelf_build);
            pause(0.001);
            [abt, abth, abtim] = neuroelf('ne_about', 0, 0, true);
            abt.CloseRequestFcn = '';
            abt.Color = [1, 1, 1];
            abt.Name = ['NeuroElf v' neuroelf_version];
            abt.Position(4) = 230;
            abth.LB_Progress.Visible = 'off';
            abth.IM_Splash.Position(2) = 4;
            feval(get(abtim, 'StopFcn'), 0, 0, true);
            abt.Visible = 'on';
            drawnow;
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            abt = struct('Delete', []);
        end

        % make sure neuroelf library is available (ne_methods is populated)
        n = neuroelf;

        % open main figure
        fprintf('figure...');
        pause(0.001);
        fMainFig = xfigure([neuroelf_path('tfg') '/neuroelf.tfg']);

        % and initialize UI/config
        n.ne_initwindow(fMainFig);

        % remove splash
        try
            abt.Delete;
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
        end

    % handle errors
    catch ne_eo;
        error('neuroelf:GUI:noGUI', 'Error creating main figure: %s.', ne_eo.message);
    end
end

% check if already in callback and don't allow additional calls then
if ne_gcfg.c.incb
    return;
end

% external and UI calls (for controls that don't allow @ne_XXX Callbacks)
if nargin > 0 && ischar(varargin{1}) && ~isempty(varargin{1}) && ...
    isfield(ne_gcfg.c.callbacks, lower(varargin{1}(:)'))

    % pass on call
    try

        % with output arguments
        if nargout > 0
            [varargout{1:nargout}] = ...
                feval(ne_gcfg.c.callbacks.(lower(varargin{1}(:)')), 0, 0, varargin{2:end});

        % without outputs
        else
            feval(ne_gcfg.c.callbacks.(lower(varargin{1}(:)')), 0, 0, varargin{2:end});
        end

    % handle errors here
    catch ne_eo;
        rethrow(ne_eo);
    end

% no arguments
elseif nargin == 0

    % just main UI page
    ne_methods.ne_showpage(0, 0, 1);
end
