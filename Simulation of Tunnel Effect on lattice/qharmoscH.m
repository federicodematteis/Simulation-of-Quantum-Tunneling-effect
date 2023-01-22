function [H,x,Htri] = qharmoscH(N,BC)
% [H,x,Htri] = dqharmosc(nx,BC) returns the "best" position-representation
% Hamiltonian matrix H and the corresponding 1D grid x of the quantum 
% harmonic oscillator with boundary conditions BC, which
% can be PBC (periodic, default), DBC (Dirichlet) or NBC (Neumann). 
% It also returns for comparison the tridiagonal approximation Htri.
% If nargout==0 a figure is opened with a plot of the two spectra.
% Example calls:  qharmoscH;
%                 qharmoscH(600);
%                 H = qharmoscH(1024,'DBC');

if nargin < 2 || isempty(BC)
    BC = 'PBC';
end
if nargin < 1 || isempty(N)
    N = 1024;
end

maxV = 0.5*(N/2)^2;  % not exact but very close
maxT = pi^2/2;       % not exact (for DBC and NBC) but very close
a = (maxT/maxV)^(1/4); % "optimal" lattice spacing

%funciton handle per definire il potenziale
V = @(x)0.5*x.^2;

% H = hamiltonian1D(N,a,V,BC) returns the Hamiltonian matrix H 
% with full infinite-order Laplacian and potential V on a uniform 1D grid 
% with N points and lattice spacing a. 
%N = number of points in the grid
[H,x,Htri] = hamiltonian1D(N,a,V,BC);

if nargout == 0
    figure
    E = eig(H);
    Etri = eig(Htri);
    plot(E)
    hold on
    plot(Etri)
    legend('H spectrum','Htri spectrum','Location','northwest')
    hold off
end