function figObject = ezSPM_Figure(FigureName, varargin)
% function ezSPM_Figure(FigureName, LinkTag, LinkObject, ...)
%
% try to open a figure named by FigureName, if FieldLinkCont is
% given it must be either an ICCRIniFile object handle or a
% struct with only ICCRIniFile handles for values
%
% FORMAT:       figObject = ezSPM_Figure(TFGFilename, IniObjectTag, IniObject)
%
% Input fields:
%
%       TFGFilename   file name of TFG file, ('ezSPM_' will be prepended)
%       IniObjectTag  tag that identifies the ini object (i.e. the first
%                     cell in any fieldlink specification line)
%       IniObject     object that will be linked to the ICCRFigure object
%
% any number of objects can be linked like this:
%
% fobj = ezSPM_Figure(TFGFile, 'tag1', obj1, 'tag2', obj2, ...)
% 
%-----------------------------------------------------------------------
% ICCR CNS - Interdisciplinary Center for Clinical Research
%            Central Nerve System
% author: Jochen Weber, June 2005, v2.0.1
% build:  5063009
%-----------------------------------------------------------------------
% ezSPM_Figure.m

% argument check
if nargin < 1 | ~ischar(FigureName)
    error( ...
        'ezSPM:InvalidArgument', ...
        'ezSPM_Figure requires the first argument to be of type CHAR.' ...
    );
end

% all ezSPM functions share this one global variable (struct)
global ezSPM_glob;

% issue ezSPM_GlobCheck to ensure the settings are loaded
ezSPM_GlobCheck;

% is FigureName no filename ?
if ~any(FigureName==filesep)
    
    % do we have the 'ezSPM_' heading ?
    if isempty(strfind(FigureName,'ezSPM'))
        FigureName = ['ezSPM_' FigureName];
    end
    
    FigureName = [ezSPM_glob.path.figures filesep FigureName];
    
    % do we have the '.tfg' ending ?
    if isempty(strfind(FigureName,'.tfg'))
        FigureName = [FigureName '.tfg'];
    end
    
end

% if FieldLinkCont is CHAR: 'filename' then return FigureName!
if nargin > 1 & ischar(varargin{1}) & strcmp(lower(varargin{1}),'filename')
    figObject = FigureName;
    return;
end

% with or without FieldLinkCont
if nargin < 3
    
    try
        % create figure object handle
        figObject = ICCRFigure(FigureName);
        
        % test for validity
        if ~isICCRFigure(figObject,1)
            error('INVALIDFIG');
        end
        
    % handling errors
    catch
        error( ...
            'ezSPM:FigureNotLoaded', ...
            'The requested figure couldn''t be loaded: %s', ...
            FigureName ...
        );
    end
    
else
    
    % build the FieldLinkCont argument
    flc = struct;
    for ac = 1:2:(nargin-1)
        try
	    flc.(varargin{ac}) = varargin{ac+1};
	end
    end
    if length(fieldnames(flc)) > 0
        LCf = fieldnames(flc);
        for fc = 1:length(LCf)
            if ischar(flc.(LCf{fc}))
                flc.(LCf{fc}) = ICCRIniFile(flc.(LCf{fc}), 'convert');
            elseif ~isICCRIniFile(flc.(LCf{fc}),1)
                flc = rmfield(flc,LCf{fc});
            end
        end
    end
    if length(fieldnames(flc)) < 1
        error( ...
            'ezSPM:InvalidArgument', ...
            'ezSPM_Figure couldn''t find any valid FLN objects.' ...
        );
    end
    
    try
        % create figure object handle
        FieldLinkCont = struct('FieldLinkCont', flc);
        figObject     = ICCRFigure(FigureName,  FieldLinkCont);
        
        % test for validity
        if ~isICCRFigure(figObject,1)
            error('INVALIDFIG');
        end
        
    catch
        error( ...
            'ezSPM:FigureNotLoaded', ...
            'The requested figure couldn''t be loaded or linked to fields: %s', ...
            FigureName ...
        );
    end
    
end

