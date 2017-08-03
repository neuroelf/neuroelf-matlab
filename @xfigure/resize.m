function resize(xo, varargin)

% for figures
if numel(xo) ~= 1 || xo.T ~= 1
    error('neuroelf:xfigure:invalidObjectType', 'Resize is only available for figures.');
end

% build resize tree
rszuics = xo.X.figprops.rszuics;
if isempty(rszuics)
    return;
end
rsztree = struct;
for cc = length(rszuics):-1:1
    uobj = rszuics{cc};

    % remove no longer existant children
    if ~isvalidxfigure(uobj)
        xo.X.figprops.rszuics(cc) = [];
        continue;
    end

    % get resize spec
    iprop = uobj.X.loadprops;
    rspec = iprop.ResizeSpec;
    utag  = iprop.Tag;

    % try to lookup uicontrol
    try
        robj = xfigure(rspec{1});
        if ~isxfigure(robj, 1) || ~any([1, 2] == robj.T)
            error('INVALID_RESIZE_SPEC');
        end
    catch xfigerror
        neuroelf_lasterr(xfigerror);
        continue;
    end

    % set this as reference
    rsztree.(utag) = [rspec(1); {robj}; {uhnd}];
end

% check tree for self references, build a relation tree and an ordered list
rsdtree = struct;
telems = fieldnames(rsztree);
while ~isempty(telems)
    telem = telems{1};
    telems(1) = [];
    chain = {telem};
    celem = rsztree.(telem){1};
    while isfield(rsztree, celem)
        if any(strcmp(chain, celem))
            error('neuroelf:xfigure:resizeChain', ...
                'Invalid ResizeSpec given. No chaining may occur.' ...
            );
        end
        chain(end + 1) = {celem};
        telems(strcmp(telems, celem)) = [];
        celem = rsztree.(celem){1};
    end
    chain(end + 1) = {celem};
    gchain = gluetostring(chain(end:-1:1), '.');
    try
        eval(['rsdtree.' gchain '=struct;']);
    catch ne_eo;
        error( ...
            'xfigure:IllegalResizeChain', ...
            'Illegal tag(s) in resize chain found (%s).', ...
            ne_eo.message ...
        );
    end
end
chain = treeorder(rsdtree);
if ~isempty(xo.X.loadprops.Tag)
    rsztree.(xo.X.loadprops.Tag) = {};
end

% get new size until it is fixed
tsize = get(xo.H, 'Position');
pause(0.1);
nsize = get(xo.H, 'Position');
while ~all(tsize == nsize)
    pause(0.05);
    tsize = nsize;
    nsize = get(xo.H, 'Position');
end
if ~isempty(xo.X.loadprops.MinSize)
    nsize(3:4) = max(nsize(3:4), xo.X.loadprops.MinSize);
    if ~all(tsize == nsize)
        if nsize(4) > tsize(4)
            nsize(2) = nsize(2) - (nsize(4) - tsize(4));
        end
        set(xo.H, 'Position', nsize);
    end
end

% now the position is fixed and we can start resizing
for rc = 1:length(chain)
    rszelem = rsztree.(chain{rc});
    if length(rszelem) ~= 4
        continue;
    end
    rszrpos = rszelem{2};
    rszrobj = rszelem{3};
    rszrtyp = rszrobj.T;
    rsztobj = rszelem{4};
    switch (rszrtyp)
        case 1
            rsize = nsize;
            rpoint = rszrpos(1:2) .* rsize(3:4) + rszrpos(5:6);
        case 2
            try
                rsize = getprop(rszrobj, 'Position');
            catch xfigerror
                neuroelf_lasterr(xfigerror);
                continue;
            end
            rpoint = rsize(1:2) + rszrpos(1:2) .* rsize(3:4) + rszrpos(5:6);
        otherwise
            warning( ...
                'xfigure:BadRefObjectType', ...
                'Controls can only be dependent on figures or controls.' ...
            );
            continue;
    end

    tgtpos = [rpoint, 0, 0];
    for pc = [3, 4]
        xpc = pc + 4;
        switch (rszrpos(pc)), case {0}
            tgtpos(pc) = rszrpos(xpc);
        case {1}
            tgtpos(pc) = rszrpos(xpc) + rsize(pc);
        case {2}
            tgtpos(pc) = rszrpos(xpc) * rsize(pc);
        case {3}
            tgtpos(pc) = rszrpos(xpc) * nsize(pc);
        otherwise
            warning( ...
                'xfigure:BadRefPosType', ...
                'Bad positioning type: %g.', ...
                rszrpos(pc) ...
            );
            continue;
        end
    end
    try
        setprop(rsztobj, 'Position', tgtpos);
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
    end
end
