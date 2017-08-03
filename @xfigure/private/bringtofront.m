function xo = bringtofront(xo)

if xo.T > 0
    try
        redrawfig(xfigmlup(ihFPos));
        figure(xfigmlup(ihFPos));
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
    end
else
    error( ...
        'xfigure:InvalidObjectType', ...
        'BringToFront is not valid for the ROOT object.' ...
    );
end


