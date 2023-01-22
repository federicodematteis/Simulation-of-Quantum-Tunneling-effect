function [H,x,Htri] = hamiltonian1D(N,a,V,BC)
% H = hamiltonian1D(N,a,V,BC) returns the Hamiltonian matrix H 
% with full infinite-order Laplacian and potential V on a uniform 1D grid 
% with N points and lattice spacing a. 
% The boundary conditions are passed in BC, which
% can be PBC (periodic, default), DBC (Dirichlet) or NBC (Neumann). 
% [H,x] = hamiltonian1D(_) also return the x grid.
% [H,x,Htri] = hamiltonian1D(_) also return the x grid and the tridiagonal
% Hamiltonian matrix Htri. 

narginchk(2,4)
if nargin < 4, BC = 'PBC'; end
if nargin < 3 || isempty(V), V = 0; end
if N < 2
    H = []; x = []; Htri = [];
    return;
end

% the grid 
x = -(N-1)/2:(N-1)/2;
x = a*x';

% input check
if isnumeric(V)
    if isscalar(V)
        V = V*ones(size(x));
    end
    if length(V) ~= N
        error('numeric potential has wrong size')
    end
elseif isa(V,'function_handle')
    V = V(x);
else
    error('potential has wrong type')
end

D = diag(ones(N-1,1),1);

switch BC
    case 'DBC'
        k = (pi/(N+1))*(1:N)';
        %L = sinft(diag(k.^2)*sinft(eye(N)));
        L = idst(diag(k.^2)*dst(eye(N)));
    case 'NBC'
        k = (pi/N)*(0:N-1)';
        %L = icosft(diag(k.^2)*cosft(eye(N)));
        L = idct(diag(k.^2)*dct(eye(N)));
        D(1,1) = 1/2; D(N,N) = 1/2;
    case 'PBC'
        n = floor(N/2);
        nn = floor((N-1)/2);
        k = (2*pi/N)*(-n:nn)';
        %L = isfft(diag(k.^2)*sfft(eye(N)));
        L = ifft(diag(fftshift(k.^2))*fft(eye(N)));
        D(N,1) = 1;
    otherwise
        error('unknown type of boundary conditions')
end

L = real(L+L')/2;
H = L/2/a^2 + diag(V);

if nargout == 3
    Htri = (eye(N) - (D+D')/2)/a^2 + diag(V);   
end
