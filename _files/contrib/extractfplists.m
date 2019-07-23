% script to extract faculty pairings from text files and write out
% several CSV files (for Excel usage)

% list of text files
institutes = { ...
    'CalTech', ...
    'Columbia', ...
    'Harvard', ...
    'MIT', ...
    'Princeton', ...
    'Stanford', ...
    'UCL', ...
    'UCSF'};

% iterate over text files
for fc = 1:numel(institutes)
    
    % grab information
    [numids, authors, singleids, pairids, fig, ax] = ...
        facultypairings([institutes{fc} '.txt']);
    
    % save image
    print(fig, [institutes{fc} '.png'], '-dpng', '-r300');
    
    % (1) turn this into a CSV with numbers of (pairwise) publications
    pwpubs = cell(numel(authors) + 1, numel(authors) + 1);
    pwpubs{1} = 'Authors';
    pwpubs(2:end, 1) = authors(:);
    pwpubs(1, 2:end) = authors(:)';
    for a1 = 1:numel(authors)
        pwpubs{a1+1, a1+1} = 'NA';
        for a2 = 1:(a1 - 1)
            pwpubs{a1+1, a2+1} = sprintf('%d', pairids(a1, a2));
            pwpubs{a2+1, a1+1} = pwpubs{a1+1, a2+1};
        end
    end
    pwpubss = sprintf([repmat('%s,', 1, numel(authors)) '%s\n'], pwpubs{:});
    fid = fopen([institutes{fc} '_matrix.csv'], 'w');
    fwrite(fid, uint8(pwpubss), '*uint8');
    fclose(fid);
    
    % (2) write out the indiviual authors publication IDs as a list
    % (a) for collaborative pairs
    for a1 = 1:numel(authors)
        pwpubs{a1+1, a1+1} = '';
        for a2 = 1:(a1 - 1)
            pwpubs{a1+1, a2+1} = sprintf('%d ', intersect(singleids{a1}, singleids{a2}));
            pwpubs{a2+1, a1+1} = pwpubs{a1+1, a2+1};
        end
    end
    pwpubss = sprintf([repmat('%s,', 1, numel(authors)) '%s\n'], pwpubs{:});
    fid = fopen([institutes{fc} '_matrix_ids.csv'], 'w');
    fwrite(fid, uint8(pwpubss), '*uint8');
    fclose(fid);
    
    % (b) for individual authors
    spubs = cell(numel(authors) + 1, 2);
    spubs(1, :) = {'Authors', 'PubMed IDs'};
    spubs(2:end, 1) = authors(:);
    for a1 = 1:numel(authors)
        spubs{a1+1, 2} = sprintf('%d,', singleids{a1});
        spubs{a1+1, 2}(end) = [];
    end
    spubs = spubs';
    spubss = sprintf('%s,%s\n', spubs{:});
    fid = fopen([institutes{fc} '_author_ids.csv'], 'w');
    fwrite(fid, uint8(spubss), '*uint8');
    fclose(fid);
    
    % delete figure
    delete(fig);
end
