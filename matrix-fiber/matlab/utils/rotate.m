a = [1 0 0
     0 0 0
 0 0 0] ;

    
 
nx = 0;
ny = 1;
nz = 0;
nx = nx/sqrt(nx.^2+ny.^2+nz.^2)
ny = ny/sqrt(nx.^2+ny.^2+nz.^2)
nz = nz/sqrt(nx.^2+ny.^2+nz.^2)
rotate_matrix = [nx,-ny,nx*nz;
    ny nx -ny*nz;
    nz 0 nx^2+ny^2];

b =  rotate_matrix'*a*rotate_matrix
% c =  rotate_matrix*b*rotate_matrix'
% c-a



