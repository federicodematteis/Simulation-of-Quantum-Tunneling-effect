function VquarticH

N=1e3;
V = @(x)x.^4/24;

% in the position representation
maxkinE = pi^2/2;
optimdx = (maxkinE/V(N/2))^(1/6);
H = hamiltonian1D(N,optimdx,V,'DBC');
e = eig(H);

% in the number representation
a = diag(sqrt(1:N-1),1);
q = (a+a')/sqrt(2);
H0 = diag((0:N-1)+1/2); 
x = eig(q);
X=max(x);
p=(a-a')*1i/sqrt(2);
momentum=eig(p);

%, (w*max(p))^2/2 = (max(x)/w)^4/24, max(p)=max(x)
w = (max(x)^2/12)^(1/6);
H1 = w^2*(H0-q^2/2)+(q/w)^4/24;
e(:,2) = eig(H1);