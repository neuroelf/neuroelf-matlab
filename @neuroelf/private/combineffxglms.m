function cffx = combineffxglms(glms)

% load GLMs
glms = glms(:);
lglms = false(size(glms));
for c = 1:numel(glms)
    if ~isxff(glms{c}, 'glm')
        try
            glms{c} = xff(glms{c});
            lglms(c) = true;
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            warning( ...
                'neuroelf:xfferror', ...
                'Error loading GLM file: %s.', ...
                ne_eo.message ...
            );
            clearxffobjects(glms(lglms));
            return;
        end
    end
end

% until everything is combined
while numel(glms) > 1

    % combine in pairs
    for c = 1:2:numel(glms)

        % attempt FFX-Join
        try
            if numel(glms) > c
                jglm = glms{c}.JoinFFX(glms{c+1});
                if lglms(c)
                    glms{c}.ClearObject;
                end
                if lglms(c+1)
                    glms{c+1}.ClearObject;
                end
                glms{c} = jglm;
                lglms(c) = true;
                lglms(c+1) = false;
            end
        catch ne_eo;
            neuroelf_lasterr(ne_eo);
            warning( ...
                'neuroelf:xfferror', ...
                'Error FFX-joining GLMs: %s.', ...
                ne_eo.message ...
            );
            clearxffobjects(glms(lglms));
            return;
        end
    end
    cglms = false(size(glms));
    cglms(2:2:end) = true;
    clearxffobjects(glms(lglms & cglms));
    glms(2:2:end) = [];
    lglms = true(size(glms));
end

% return object
cffx = glms{1};
