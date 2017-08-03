function os = makeostruct(itype)
os = cell2struct(cell(1, 1, 8), {'callbacks', 'deletefcn', 'figtype', ...
    'figprops', 'loadprops', 'prevprops', 'timeclick', 'uicprops'}, 3);
os.callbacks = cell(1, 4);
os.prevprops = struct;
if itype == 1
    os.figprops = cell2struct(cell(1, 1, 12), {'cpage', 'egroups', ...
        'lgroups', 'lilookup', 'linkcont', 'linkspec', 'llookup', ...
        'pages', 'rgroups', 'rszuics', 'sgroups', 'vgroups'}, 3);
end
