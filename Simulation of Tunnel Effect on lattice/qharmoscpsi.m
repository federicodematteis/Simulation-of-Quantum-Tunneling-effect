function [psi,x,w] = qharmoscpsi(n,x)
% psi = qharmoscpsi(n) returns the first n eigenfunctions of the quantum
% hamonic oscillator H=p^2/2+q^2/2, evaluated over the zeros of H_n, the
% Hermite polynomial of order n, using the Golub-Welsch algorithm.
% psi is a matrix whose columns are the wavefunctions values.
%E
% [psi,x,w] = qharmoscpsi(n) returns in column array x the Gauss-Hermite 
% nodes of order n, that is the zeros of H_n and, in column array w, the 
% Gauss integration weights for the measure exp(-y^2)*dy multiplied by
% exp(x.^2) (the so-called scaled weights).
%
% psi = qharmoscpsi(n,x) returns the first n eigenfunctions evaluated in 
% arbitrary 1D array x using the recursion rule of Hermite polynomials. 
% For instance, psi(:,1)=exp(-x(:).^2/2)/pi^(1/4) is the wavefunction 
% of the ground state.
% If x is the set of zeros of H_n, then [psi,~,w] = qharmoscpsi(n,x)
% returns also (much) more accurate values of w for larger x. From these
% values, the Gauss integration weights for the measure exp(-y^2)*dy
% can be accurately computed by dividing by exp(x.^2).
% The case x = [] is treated as the zeros of H_n.
%
% N.B.: double precision arithmetics has troubles when exp(x^2) gets close
% to realmax, which happens for x near 26.64, but with the recursion method
% w can be computed with 1e-8 accuracy up to x=38.169918....when n = 765.

narginchk(1,2)
if nargin == 1
    a = diag(sqrt(1:n-1),1); % truncated annihilation operator
    q = (a+a')/sqrt(2);      % truncated position operator (companion matrix)
    [u,x] = eig(q);
    x = diag(x);
    x = (x - flip(x))/2; % force odd symmetry;
    w = u(1,:).^2; %Gauss weights for the measure exp(-y^2)*dy/sqrt(pi)
    u = u./sign(u(1,:)); % this fixes the sign at different values of x
    % we need to force even/odd symmetry;
    u(1:2:end,:) = (u(1:2:end,:) + fliplr(u(1:2:end,:)))/2;
    u(2:2:end,:) = (u(2:2:end,:) - fliplr(u(2:2:end,:)))/2;
    % u./sqrt(w) is the matrix of Hermite polynomials evaluated at the
    % zeros of H_n (along the rows), but for large n it would work only 
    % for small enough x due to the finite precision. 
    % If we compute the wavefunctions the problem is less serious but is
    % still there for large enough x and large enough orders.
    w = (w + flip(w))/2; % force even symmetry;
    w = sqrt(pi)*exp(x.^2 + log(w'));
    psi = (u./sqrt(w'))';
    return
end

if isrow(x), x = x'; end
if isempty(x)
    a = diag(sqrt(1:n-1),1);
    q = (a+a')/sqrt(2);
    x = eig(q);
end

psi0 = exp(-x.^2/2)/pi^(1/4);
psi = psi0;
if n == 0, return, end

psi1 = sqrt(2)*x.*psi0;
psi = [psi0,psi1];
if n == 1, return, end

psi = [psi,zeros(length(x),n-2)];
for j = 2:n-1
    psi(:,j+1) = sqrt(2/j)*(x.*psi1 - sqrt((j-1)/2)*psi0);
    psi0 = psi1;
    psi1 = psi(:,j+1);
end

w = [];
% This is an alternative, more accurate way to compute w.
% Still, in double precision, NaN appear at large enough x when n > 765
if n == length(x)
    psi_n = sqrt(2/n)*(x.*psi1 - sqrt((n-1)/2)*psi0);
    if norm(psi_n,inf) < 1e-10  % input x are accurate zeros of H_n
        u = psi./sqrt(sum(psi.^2,2));
        w = sqrt(pi)*exp(2*log(u(:,1))+x.^2);
    end
end
if isempty(w)
    warning('cannot compute weigths with the given inputs')
end

    