N = 2.24e3;
a = 0.1;
%changed in the last lesson
%a = 0.1;

tau = 0.01;
%nt = 6e3;
%nt = 8e3;
nt = 1e4;
x = -(N-1)/2:(N-1)/2; 
x = a*x';

%%
shg; clf
set(gcf,'position',[900,600,800,400]);
%V = @(x)exp(-x.^2/8);
V = @(x)0.6*double(x.^2 < 6);
hV = plot(x,V(x),'linewidth',2);  %plot del potenziale 
hold on
xlim([x(1) x(end)]);
axis(axis)

%%
% k0 = 1.5;
%qui gli dico quale momento usare 
k0 = 1;
sigma = 20;
psi0 = exp(1i*k0*x).*exp(-(x+50).^2/sigma^2/2);
psi0 = psi0/norm(psi0);
psi = psi0;   
s = 150;
%serve a non riplottare se esiste gia il graphic handle
%che corrisponde al disegno della funzione d'onda.
if exist('hf','var') && isgraphics(hf)
    hf.YData = s*abs(psi).^2;        %plot dello stato sul grafico
else
    hf = plot(x,s*abs(psi).^2,'linewidth',2);
end
ht = title(sprintf('t = %-5.2f',0));

%definisce k della grid
k = (pi/(N+1))*(1:N)';
%definisce l'energia cinetica
kinE = (k/a).^2/2; %diagonale nello spazio di Fourier

shg
disp('press any key to start execution'); pause

for j=1:nt
    psi = exp(-1i*V(x)*tau/2).*psi;
    psi = sinft(exp(-1i*kinE*tau).*sinft(psi));
    psi = exp(-1i*V(x)*tau/2).*psi;
    if mod(j,10) == 0
        hf.YData = s*abs(psi).^2;
        ht.String = sprintf('t = %-5.2f',j*tau);
        drawnow;
    end
end

[H,~,Htri] = hamiltonian1D(N,a,V,'DBC');

return

%%
psi = psi0;
shg
U = expm(-1i*H*tau);
disp('press any key to start execution'); pause

for j=1:nt
    psi = U*psi; 
    if mod(j,10) == 0
        hf.YData = s*abs(psi).^2;
        ht.String = sprintf('t = %-5.2f',j*tau);
        drawnow;
    end
end

%%
n = 30;
yV = max(V(x))/2
psi = (u(:,n)-1i*u(:,n+1))/sqrt(2);
hf.YData = sqrt(s)*real(psi)/2 + max(V(x))/2;
hf2 = plot(x,sqrt(s)*imag(psi)/2 + yV,'linewidth',2);
%hf.YData = s*abs(psi).^2;

shg
disp('press any key to start execution'); pause

for j=1:nt
    psi = exp(-1i*e(n)*tau)*psi; 
    if mod(j,10) == 0
        hf.YData = sqrt(s)*real(psi)/2 + max(V(x))/2;
        hf2.YData = sqrt(s)*imag(psi)/2 + max(V(x))/2;
        %hf.YData = s*abs(psi).^2;
        ht.String = sprintf('t = %-5.2f',j*tau);
        drawnow;
        
    end
end

    