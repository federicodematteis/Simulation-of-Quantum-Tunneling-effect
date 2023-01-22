function L = laplacian1D(N,p,BC)
% L = laplacian1D(N,p,BC) returns the discretized sparse 1D laplacian 
% of order 2*p+2 and boundary conditions BC, that can be:
% 'PBC'   Periodic BC, the default;
% 'DBC0'  all zeros Dirichlet BC) namely f(j)=0 for j<=0 and j>=N+1;
% 'DBC'   Dirichlet BC: f(-j)=-f(j) and f(N+1+j)=-f(N+1-j) for j>=0;
% 'DSTI'  DBC through discrete sine transform of type I;
% 'NBC'   Neumann BC: f(1-j)=f(j) and f(N+j)=f(N-j+1) for j>0;
% 'DCTII' NBC through discrete cosine transform of type II,
% If p = inf, the dense Laplacian with exact k^2 spectrum is returned. In
% this case 'DBC0' is treated as 'DBC'. 

