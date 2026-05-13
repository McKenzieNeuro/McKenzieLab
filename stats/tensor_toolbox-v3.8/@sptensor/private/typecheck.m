function typecheck(type)

if ~strcmp(type,'sparse') && ~strcmp(type,'incomplete')
    error('Type must be ''sparse'' or ''incomplete''.')
end