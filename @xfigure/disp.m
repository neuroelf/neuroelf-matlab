function disp(xo)
%XFIGURE::DISP  Display information for an object.

% Version:  v1.1
% Build:    16041817
% Date:     Apr-18 2016, 5:25 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2011, 2016, Jochen Weber
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

% global reference storage
global ne_methods xfigmlup xfigsngl xfigures;

% iterate over objects and report last error if any
for oc = 1:numel(xo)
    try
        % switch over type
        rqobj  = xo(oc);
        switch rqobj.T
            case 0
                lTyp = 'ROOT';
            case 1
                lTyp = 'Figure';
            case 2
                lTyp = 'UIControl';
            case 3
                lTyp = 'UIMenu';
            case 4
                lTyp = 'UIContextMenu';
            otherwise
                lTyp = 'UNKNOWN';
        end
        mentxt = '';

        % root object
        switch rqobj.T, case 0

            % prepare output string
            tgsc = '';

            % garbage collection
            xfiggc();
            
            % iterate over figures
            mltypes = lower(get(xfigmlup, 'Type'));
            if ~iscell(mltypes)
                mltypes = cellstr(mltypes);
            end
            myfigures = find(strcmp(mltypes, 'figure'));
            for fc = 1:numel(myfigures)

                % get shortcut and test availability
                mylfigure = xfigures(myfigures(fc));
                if ~isvalidxfigure(mylfigure)
                    warning('neuroelf:xfigure:figureDisappeared', ...
                        'Figure disappeared from class control.')
                    continue;
                end

                % get group names
                tt = ['  Groups of figure ''' mylfigure.X.loadprops.Title ...
                      ''' (figure no. ' num2str(fc) '):' char(10) char(10)];
                if isstruct(mylfigure.X.figprops.egroups)
                    egs = fieldnames(mylfigure.X.figprops.egroups);
                    egss = cell(1, numel(egs));
                    for gc = 1:numel(egs)
                        egc  = mylfigure.X.figprops.egroups.(egs{gc});
                        luic = numel(intersect(egc, xfigures(strcmp(mltypes, 'uicontrol'))));
                        luim = numel(intersect(egc, xfigures(strcmp(mltypes, 'uimenu'))));
                        egss{gc} = sprintf('%40s: %d UIControls, %d UIMenus%s', ...
                            egs{gc}, luic, luim, char(10));
                    end
                    ego = ['    EGroups:' char(10) char(10) ne_methods.gluetostring(egss, '')];
                else
                    ego = '';
                end
                if isstruct(mylfigure.X.figprops.rgroups)
                    rgs = fieldnames(mylfigure.X.figprops.rgroups);
                    rgos = cell(1, numel(rgs));
                    for gc = 1:numel(rgs)
                        rgc = mylfigure.X.figprops.rgroups.(rgs{gc});
                        rgos{gc} = sprintf('%40s: %d UIControls%s', rgs{gc}, numel(rgc), char(10));
                    end
                    rgo = ['    RGroups:' char(10) char(10) ne_methods.gluetostring(rgos, '')];
                else
                    rgo = '';
                end
                if isstruct(mylfigure.X.figprops.sgroups)
                    sgs = fieldnames(mylfigure.X.figprops.sgroups);
                    sgos = cell(1, numel(sgs));
                    for gc = 1:numel(sgs)
                        sgc = mylfigure.X.figprops.sgroups.(sgs{gc});
                        sgos{gc} = sprintf('%40s: %d UIControls%s', sgs{gc}, numel(sgc), char(10));
                    end
                    sgo = ['    SGroups:' char(10) char(10) ne_methods.gluetostring(sgos, '')];
                else
                    sgo = '';
                end
                if isstruct(mylfigure.X.figprops.vgroups)
                    vgs = fieldnames(mylfigure.X.figprops.vgroups);
                    vgos = cell(1, numel(vgs));
                    for gc = 1:numel(vgs)
                        vgc = mylfigure.X.figprops.vgroups.(vgs{gc});
                        vgos{gc} = sprintf('%40s: %d UIControls%s', vgs{gc}, numel(vgc), char(10));
                    end
                    vgo = ['    VGroups:' char(10) char(10) ne_methods.gluetostring(vgos, '')];
                else
                    vgo = '';
                end

                % add figure output to list
                tgsc = sprintf('%s%s%s%s%s%s',tgsc, tt, ego, rgo, sgo, vgo);
            end

            % what tags are currently in use
            uto = ['  Tags for find used in xfigure base class:' char(10) char(10)];
            uts = fieldnames(xfigsngl.tags);
            utss = cell(1, numel(uts));
            for utc = 1:numel(uts)
                utag = uts{utc};
                utt  = xfigsngl.tags.(utag);
                if ~isvalidxfigure(utt)
                    xfigsngl.tags = rmfield(xfigsngl.tags, utag);
                    continue;
                end
                uttt = utt.T + 1;
                if uttt > 0 && uttt < numel(xfigsngl.objtypel)
                    tlTyp = xfigsngl.objtypel{uttt};
                    if uttt > 1
                        utag = utag(5:end);
                    end
                else
                    tlTyp = 'xfigure_unknown_type';
                end
                utss{utc} = sprintf('%36s: tag of a subtype %s object%s', ...
                    utag, tlTyp, char(10));
            end
            uto = [uto ne_methods.gluetostring(utss, '')];

            % format output with some other properties
            Props = sprintf([ ...
                        '    Figures         active:  %4.0f\n' ...
                        '    UIControls      active:  %4.0f\n' ...
                        '    UIMenus         active:  %4.0f\n' ...
                        '    UIContextMenus  active:  %4.0f\n' ...
                        '    xfigure Tags    active:  %4.0f\n' ...
                        '\n%s%s'], ...
                    sum(strcmpi(mltypes, 'figure')), ...
                    sum(strcmpi(mltypes, 'uicontrol')), ...
                    sum(strcmpi(mltypes, 'uimenu')), ...
                    sum(strcmpi(mltypes, 'uicontextmenu')), ...
                    numel(fieldnames(xfigsngl.tags)), ...
                    tgsc, uto);

        % figures
        case 1
            lProps = xfigsngl.outfig;

        % uicontrols
        case 2
            lProps = xfigsngl.outuic;

        % uimenus
        case 3
            lProps = xfigsngl.outuim;
            mentxt = menutext(rqobj.H);

        % uicontextmenus
        case 4
            lProps = xfigsngl.outuix;
            mentxt = menutext(rqobj.H);

        % new object
        case -1
            lTyp = 'placeholder';
            Props = 'none';

        % otherwise raise an error
        otherwise
            error( ...
                'xfigure:InvalidObjectType', ...
                'Can''t display information for unknown object type.' ...
            );
        end

        % collection information for non-ROOT objects
        if rqobj.T > 0
            fns    = fieldnames(lProps);
            if xfigsngl.mlversion < 804
                Props  = sprintf([ ...
                         '%24s: %.8f\n', ...
                         '%24s: %s\n\n'], ...
                         'MATLAB GUI handle', double(rqobj.H), ...
                         'Tag (OnLoad)', rqobj.X.loadprops.Tag);
            else
                Props = sprintf([ ...
                        '%24s: object\n', ...
                        '%24s: %s\n\n'], ...
                        'Matlab GUI handle', ...
                        'Tag (OnLoad)', rqobj.X.loadprops.Tag);
            end
            for fnc = 1:numel(fns)
                pnow = [];
                for rqc = rqobj(:)'
                    try
                        pnow = getprop(rqc, lProps.(fns{fnc}));
                        break;
                    catch ne_eo;
                        neuroelf_lasterr(ne_eo);
                    end
                end
                if (~ischar(pnow) || numel(pnow) ~= numel(pnow)) && ~isempty(pnow)
                    if any(strcmpi(class(pnow), ...
                        {'double', 'int16', 'int32', 'int64', 'int8', 'logical', ...
                         'single', 'uint16', 'uint32', 'uint64', 'uint8'})) && ...
                        numel(pnow) <= 8
                        pnow = any2ascii(pnow);
                    elseif isa(pnow, 'function_handle')
                        pnow = sprintf('function handle: %s', func2str(pnow));
                    else
                        psiz = sprintf('%.0fx', size(pnow));
                        pcls = class(pnow);
                        if isstruct(pnow)
                            pfld = sprintf(' with %.0f field(s)', numel(fieldnames(pnow)));
                        else
                            pfld = '';
                        end
                        pnow = sprintf('<%s %s%s>', psiz(1:end-1), pcls, pfld);
                    end
                else
                    pnow = pnow(:)';
                end
                if ~isempty(pnow)
                    Props = sprintf('%s%24s: %s\n', Props, fns{fnc}, pnow);
                end
            end

            if ~isempty(mentxt)
                Props = [Props char(10) 'Menu tree:' char(10) mentxt];
            end
        end

        % display output to console
        disp(['xfigure object of type ' lTyp '. Properties:' char(10) Props]);
    catch ne_eo;
        rethrow(ne_eo);
    end
end
