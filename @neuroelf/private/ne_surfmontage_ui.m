% FUNCTION ne_surfmontage_ui: update UI based on callback
function ne_surfmontage_ui(varargin)

% Version:  v1.1
% Build:    16053009
% Date:     May-30 2016, 9:46 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2015, 2016, Jochen Weber
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

% global variable
global ne_gcfg;

% get tags
sm = ne_gcfg.h.SurfMontage;
cc = ne_gcfg.fcfg.SurfMontage;
hFig = sm.SMFig;
tags = sm.h;

% current config
fcc = fieldnames(cc.cc);
cci = 0;
if ~isempty(fcc)
    cci = tags.DD_surfmontage_configs.Value;
    ccc = cc.cc.(fcc{cci});
else
    ccc = {};
end

% depending on what to do
switch (lower(varargin{3}))

    % add configuration
    case 'addconfig'

        % get configuration name
        cfgname = inputdlg({'Please enter the name for the new configuration:'}, ...
            'NeuroElf - user input', 1, {'surface montage'});
        if isempty(cfgname) || ~iscell(cfgname) || ~ischar(cfgname{1}) || ...
            isempty(ddeblank(cfgname{1}))
            return;
        end
        cfgname = ddeblank(cfgname{1});

        % add configuration

    % background color selection
    case 'bgcolor'
        if cci == 0
            return;
        end

        % request new color
        newcol = colorpicker(ccc{3}, {'Surface montage background color'});
        if ~isequal(newcol, ccc{3})
            ne_gcfg.fcfg.SurfMontage.cc.(fcc{cci}){3} = newcol;
            tags.BT_surfmontage_bgcol.CData = ...
                repmat(reshape(uint8(newcol), [1, 1, 3]), [24, 36, 1]);
        end

    % check field contents and update config
    case 'checkfield'
        if nargin < 4 || ~ischar(varargin{4}) || isempty(varargin{4})
            return;
        end

        % which field to check
        switch (lower(varargin{4}(:)'))

            % config wide fields
            case {'imagex', 'imagey'}
                if cci == 0 || isempty(ccc)
                    return;
                end

                % get values
                w = str2double(tags.ED_surfmontage_imagex.String);
                h = str2double(tags.ED_surfmontage_imagey.String);
                if numel(w) ~= 1 || isinf(w) || isnan(w) || w < 16 || w > 10240 || ...
                    numel(h) ~= 1 || isinf(h) || isnan(h) || h < 16 || h > 10240
                    w = ccc{2}(1);
                    h = ccc{2}(2);
                else
                    w = round(w);
                    h = round(h);
                end
                tags.ED_surfmontage_imagex.String = sprintf('%d', w);
                tags.ED_surfmontage_imagey.String = sprintf('%d', h);
                ne_gcfg.fcfg.SurfMontage.cc.(fcc{cci}){2} = [w, h];

            % element wide fields
            case {'azimuth', 'elemx', 'elemy', 'infforce', 'infniter', ...
                  'offsetx', 'offsety', 'smforce', 'smniter', ...
                  'time', 'transx', 'transy', 'zenith', 'zoom'}
                if cci == 0 || isempty(ccc) || isempty(ccc{5})
                    return;
                end
                ecf = fieldnames(cc.ec);
                eci = tags.DD_surfmontage_element.Value;
                ecf = ecf{eci};
                elm = cc.ec.(ecf);
                if numel(elm) ~= 15
                    return;
                end
                val = str2double(tags.ED_surfmontage_smniter.String);
                if numel(val) ~= 1 || isinf(val) || isnan(val) || val < 0
                    val = elm{2};
                else
                    val = round(min(val, 5000));
                end
                elm{2} = val;
                val = str2double(tags.ED_surfmontage_smforce.String);
                if numel(val) ~= 1 || isinf(val) || isnan(val) || val <= 0 || val >= 1
                    val = elm{3};
                else
                    val = 0.001 * round(1000 * val);
                end
                elm{3} = val;
                val = str2double(tags.ED_surfmontage_infniter.String);
                if numel(val) ~= 1 || isinf(val) || isnan(val) || val < 0
                    val = elm{4};
                else
                    val = round(min(val, 1000));
                end
                elm{4} = val;
                val = str2double(tags.ED_surfmontage_infforce.String);
                if numel(val) ~= 1 || isinf(val) || isnan(val) || val <= 0 || val >= 1
                    val = elm{5};
                else
                    val = 0.001 * round(1000 * val);
                end
                elm{5} = val;
                val = str2double(tags.ED_surfmontage_transx.String);
                if numel(val) ~= 1 || isinf(val) || isnan(val) || val <= -256 || val >= 256
                    val = elm{6};
                end
                elm{6} = round(val);
                val = str2double(tags.ED_surfmontage_transy.String);
                if numel(val) ~= 1 || isinf(val) || isnan(val) || val <= -256 || val >= 256
                    val = elm{7};
                end
                elm{7} = round(val);
                val = str2double(tags.ED_surfmontage_azimuth.String);
                if numel(val) ~= 1 || isinf(val) || isnan(val)
                    val = elm{8};
                else
                    val = mod(val, 360);
                end
                elm{8} = round(val);
                val = str2double(tags.ED_surfmontage_zenith.String);
                if numel(val) ~= 1 || isinf(val) || isnan(val) || val < -90 || val > 90
                    val = elm{9};
                end
                elm{9} = round(val);
                val = str2double(tags.ED_surfmontage_zoom.String);
                if numel(val) ~= 1 || isinf(val) || isnan(val) || val < 0.2 || val > 5
                    val = elm{10};
                else
                    val = 0.001 * round(1000 * val);
                end
                elm{10} = val;
                val = str2double(tags.ED_surfmontage_time.String);
                if numel(val) ~= 1 || isinf(val) || isnan(val) || val < 0 || val > 30
                    val = elm{11};
                else
                    val = 0.001 * round(1000 * val);
                end
                elm{11} = val;
                val = str2double(tags.ED_surfmontage_elemx.String);
                if numel(val) ~= 1 || isinf(val) || isnan(val) || ...
                    val < 0 || val > 2560 || (val >= 1 && val < 16)
                    val = elm{12};
                end
                if val > 1
                    elm{12} = round(val);
                else
                    elm{12} = val;
                end
                val = str2double(tags.ED_surfmontage_elemy.String);
                if numel(val) ~= 1 || isinf(val) || isnan(val) || ...
                    val < 0 || val > 2560 || (val >= 1 && val < 16)
                    val = elm{13};
                end
                if val > 1
                    elm{13} = round(val);
                else
                    elm{13} = val;
                end
                val = str2double(tags.ED_surfmontage_offsetx.String);
                if numel(val) ~= 1 || isinf(val) || isnan(val) || val < 0 || ...
                   (elm{12} >= 1 && (val ~= fix(val) || val > (ccc{2}(1) - elm{12}))) || ...
                   (elm{12} < 1 && val > (1 - elm{12}))
                    if elm{12} >= 1
                        if elm{14} >= 1
                            val = elm{14};
                        else
                            val = round(elm{14} * ccc{2}(1));
                        end
                    else
                        if elm{14} < 1
                            val = elm{14};
                        else
                            val = elm{14} / ccc{2}(1);
                        end
                    end
                end
                if val >= 1
                    elm{14} = round(val);
                else
                    elm{14} = val;
                end
                val = str2double(tags.ED_surfmontage_offsety.String);
                if numel(val) ~= 1 || isinf(val) || isnan(val) || val < 0 || ...
                   (elm{13} >= 1 && (val ~= fix(val) || val > (ccc{2}(2) - elm{13}))) || ...
                   (elm{13} < 1 && val > (1 - elm{13}))
                    if elm{13} >= 1
                        if elm{15} >= 1
                            val = elm{15};
                        else
                            val = round(elm{15} * ccc{2}(2));
                        end
                    else
                        if elm{15} < 1
                            val = elm{15};
                        else
                            val = elm{15} / ccc{2}(2);
                        end
                    end
                end
                if val >= 1
                    elm{15} = round(val);
                else
                    elm{15} = val;
                end

                % store back, re-generate list (update texts), and re-set
                ne_gcfg.fcfg.SurfMontage.ec.(ecf) = elm;
                ne_surfmontage_ui(0, 0, 'genelemlist');
                ne_surfmontage_ui(0, 0, 'setconfig', eci);
                ne_surfmontage_ui(0, 0, 'setelement');

            % filename
            case 'filename'
                newfilename = tags.ED_surfmontage_filename.String;
                if isempty(newfilename) || ~ischar(newfilename) || ...
                    isempty(regexpi(ddeblank(newfilename), '^\S+\.(bmp|eps|jpg|jpeg|png|tif|tiff)$'))
                    newfilename = ne_gcfg.fcfg.SurfMontage.mc.WriteFilename;
                end
                newfilename = ddeblank(newfilename);
                ne_gcfg.fcfg.SurfMontage.mc.WriteFilename = newfilename;
                tags.ED_surfmontage_filename.String = [' ' newfilename];
                
            % copy stats bars checkbox
            case 'stbars'
                if cci == 0 || isempty(ccc)
                    return;
                end
                ne_gcfg.fcfg.SurfMontage.cc.(fcc{cci}){4} = ...
                    (tags.CB_surfmontage_stbars.Value > 0);
        end

    % generate element list
    case 'genelemlist'
        ecf = fieldnames(cc.ec);
        for eci = 1:numel(ecf)
            ecc = cc.ec.(ecf{eci});
            [epath, ecf{eci}] = fileparts(ecc{1});
            if numel(ecf{eci}) > 40
                ecf{eci} = [ecf{eci}(1:18) '...' ecf{eci}(end-17:end)];
            end
            ecf{eci} = sprintf('%s (%d/%d, %d/%d, %.3fx, %.3fs; sm=%d, inf=%d)', ecf{eci}, ...
                ecc{6:11}, ecc{2}, ecc{4});
        end
        if isempty(ecf)
            ecf = {'<no surface elements available>'};
        end
        tags.DD_surfmontage_element.String = ecf;
        
    % select element from listbox
    case 'selectelement'
        if cci == 0 || isempty(ccc) || isempty(ccc{5})
            return;
        end

        % find according piece
        ccm = ccc{5}{tags.LB_surfmontage_celems.Value};
        ecf = fieldnames(cc.ec);
        cci = findfirst(strcmpi(ecf, ccm));
        if ~isempty(cci)
            tags.DD_surfmontage_element.Value = cci;
            ne_surfmontage_ui(0, 0, 'setelement');
        end

    % set a configuration
    case 'setconfig'
        if cci == 0
            return;
        end

        % general flags
        tags.ED_surfmontage_imagex.String = sprintf('%d', ccc{2}(1));
        tags.ED_surfmontage_imagey.String = sprintf('%d', ccc{2}(2));
        tags.BT_surfmontage_bgcol.CData = ...
            repmat(reshape(uint8(ccc{3}), [1, 1, 3]), [24, 36, 1]);
        tags.CB_surfmontage_stbars.Value = double(ccc{4});

        % generate list of particle elements
        plist = ccc{5}(:);
        elist = tags.DD_surfmontage_element.String;
        edata = fieldnames(cc.ec);
        for pc = numel(plist):-1:1
            eindex = findfirst(strcmpi(plist{pc}, edata));
            if isempty(eindex)
                plist(pc) = [];
            else
                plist{pc} = elist{eindex};
            end
        end
        tags.LB_surfmontage_celems.String = plist;
        if nargin > 3 && isa(varargin{4}, 'double') && numel(varargin{4}) == 1 && ...
           ~isinf(varargin{4}) && ~isnan(varargin{4}) && varargin{4} > 0 && ...
            varargin{4} == round(varargin{4})
            tags.LB_surfmontage_celems.Value = varargin{4};
        else
            tags.LB_surfmontage_celems.Value = 1;
        end
        ne_surfmontage_ui(0, 0, 'selectelement');

    % set element dropdown -> update fields
    case 'setelement'
        ecf = fieldnames(cc.ec);
        eci = ecf{tags.DD_surfmontage_element.Value};
        elm = cc.ec.(eci);
        tags.ED_surfmontage_smniter.String = sprintf('%d', elm{2});
        tags.ED_surfmontage_smforce.String = sprintf('%g', elm{3});
        tags.ED_surfmontage_infniter.String = sprintf('%d', elm{4});
        tags.ED_surfmontage_infforce.String = sprintf('%g', elm{5});
        tags.ED_surfmontage_transx.String = sprintf('%d', elm{6});
        tags.ED_surfmontage_transy.String = sprintf('%d', elm{7});
        tags.ED_surfmontage_azimuth.String = sprintf('%d', elm{8});
        tags.ED_surfmontage_zenith.String = sprintf('%d', elm{9});
        tags.ED_surfmontage_zoom.String = sprintf('%g', elm{10});
        tags.ED_surfmontage_time.String = sprintf('%g', elm{11});
        if elm{12} == fix(elm{12})
            tags.ED_surfmontage_elemx.String = sprintf('%d', elm{12});
        else
            tags.ED_surfmontage_elemx.String = sprintf('%0.5g', elm{12});
        end
        if elm{13} == fix(elm{13})
            tags.ED_surfmontage_elemy.String = sprintf('%d', elm{13});
        else
            tags.ED_surfmontage_elemy.String = sprintf('%0.5g', elm{13});
        end
        if elm{14} == fix(elm{14})
            tags.ED_surfmontage_offsetx.String = sprintf('%d', elm{14});
        else
            tags.ED_surfmontage_offsetx.String = sprintf('%0.5g', elm{14});
        end
        if elm{15} == fix(elm{15})
            tags.ED_surfmontage_offsety.String = sprintf('%d', elm{15});
        else
            tags.ED_surfmontage_offsety.String = sprintf('%0.5g', elm{15});
        end

    % sample from VMP checkbox
    case 'smpvmp'

        % enable controls as needed
        vmpg = (tags.CB_surfmontage_samplevmp.Value > 0);
        ne_gcfg.fcfg.SurfMontage.mc.SampleVMP = vmpg;
        if vmpg
            vmpg = 'on';
        else
            vmpg = 'off';
        end
        hFig.SetGroupEnabled('SVMP', vmpg);

    % dis-/enable filename
    case 'targetfigure'
        tags.ED_surfmontage_filename.Enable = 'off';
        hFig.RadioGroupSetOne('VisMOut', 1);
    case 'targetfile'
        tags.ED_surfmontage_filename.Enable = 'on';
        hFig.RadioGroupSetOne('VisMOut', 2);
end
