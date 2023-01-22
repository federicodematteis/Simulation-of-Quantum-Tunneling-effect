function [g,a] = squarelattice(n)
% [g,a] = squarelattice(n) returns the graph object and the adjacency matrix
% of a d-dimensional square lattice with abs(n(1)) x abs(n(1))....abs(n(d))
% points. If n(j) is positive the j-directon is periodic, otherwise is
% open.

if size(n,1) > 1
    error('input must be a row vector')
end
d = length(n);

if d == 1
    pflag = n > 0;
    n = abs(ceil(n));
    a = spdiags(ones(n,1),1,n,n);
    if pflag
        a(n(1),1) = 1;
    end
    a = a + a';
    g = graph(a);
    return; 
end

[~,a] = squarelattice(n(1:end-1));
[~,b] = squarelattice(n(end));
n = abs(n);

a = kron(a,speye(n(end))) + kron(speye(prod(n(1:end-1))),b);
g = graph(a);
