classdef MeshMap < uint32
    
    % Mesh edit map identifiers.  These are used via MeshMap.MATERIAL etc.
    % The user can add more if needed.  Ideas could be to denote CMFD regions
    % (if an automated approach isn't taken).  These don't have to be used, but
    % it's handy.
    enumeration 
         MATERIAL(0)
         ASSEMBLY(1)
         PIN(2)
         FUEL(3)
         MODERATOR(4)
    end
    
end