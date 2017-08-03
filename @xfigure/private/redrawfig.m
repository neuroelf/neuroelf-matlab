function redrawfig(mlh)
%XFIGURE::PRIVATE::REDRAWFIG  Internal redraw code to force redraw.

% allow for errors
try

    % get original color
    tcol = get(mlh, 'Color');

    % set color and immediately reset color
    set(mlh, 'Color', (tcol * 0.99 + [0.002, 0.002, 0.002]), 'Color', tcol);
catch xfigerror
    neuroelf_lasterr(xfigerror);
end

% issue drawnow
drawnow;
