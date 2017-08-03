function radiogroupsetone(xo, iStr, e)
%XFIGURE::RADIOGROUPSETONE  Set one button active (enabled) in an RGroup.
%   RADIOGROUPSETONE(BUTTON) sets a single button active and all other
%   buttons in the same group inactive.

% only valid for figures and uicontrols
if numel(xo) ~= 1 || ~any([1, 2] == xo.T)
    error('neuroelf:xfigure:invalidObjectType', ...
        'RadioGroupSetOne is only valid for UIControls and figures.');
end

% for uicontrols, check Style property
if xo.T == 2 && ~strcmpi(get(xo.H, 'Style'), 'radiobutton')
    error('neureolf:xfigure:invalidObjectType', ...
        'RadioGroupSetOne is only vali for RadioButton UIControls.');
end

% get correct radio group list
if xo.T == 1
    po = xo;
else
    po = xfigure(get(xo.H, 'Parent'));
end
rgroups = po.X.figprops.rgroups;

% any radiogroups defined
if isempty(fieldnames(rgroups))
    return;
end

% for uicontrols
if xo.T == 2
    myrgroup = xo.X.loadprops.RGroup;
    if isempty(myrgroup)
        return;
    end
    if ~isfield(rgroups, myrgroup)
        error('neuroelf:xfigure:radioGroupInvalid', ...
            'The RadioButtonGroup %s is not defined in this figure.', myrgroup);
    end

    % get all controls in this group
    controls = rgroups.(myrgroup);
    for ctc = 1:numel(controls)
        try
            set(controls(ctc).H, 'Value', get(controls(ctc).H, 'Min'));
        catch xfigerror
            neuroelf_lasterr(xfigerror);
        end
    end
    try
        set(xo.H, 'Value', get(xo.H, 'Max'));
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
    end

% for figures
else
    if nargin < 3 || ~isvarname(deblank(iStr(:)')) || ~isnumeric(e) || ...
        isempty(e) || isnan(e(1)) || isinf(e(1))
        error('neuroelf:xfigure:badArgument', ...
            'For figures RadioGroupSetOne requires a group and button ID.');
    end
    dogroup = deblank(iStr(:)');

    % test group
    if ~isfield(rgroups, dogroup)
        error('neuroelf:xfigure:radioGroupInvalid', ...
            'The RadioButtonGroup %s is not defined in this figure.', dogroup);
    end
    controls = rgroups.(dogroup);
    actbutton = fix(e(1));
    if actbutton < 1 || actbutton > numel(controls)
        error('neuroelf:xfigure:badArgument', ...
            'The RadioButtonGroup doesn''t have a button #%d.', actbutton);
    end

    % iterate
    for ctc = 1:numel(controls)
        try
           if ctc ~= actbutton
               set(controls(ctc).H, 'Value', get(controls(ctc).H, 'Min'));
           else
               set(controls(ctc).H, 'Value', get(controls(ctc).H, 'Max'));
           end
        catch xfigerror
            neuroelf_lasterr(xfigerror);
        end
    end
end
