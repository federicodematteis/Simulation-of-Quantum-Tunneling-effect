function lattice
format long
n1=[-12];
n2=[12];
g=squarelattice(n1);
h=squarelattice(n2);
plot(g)
hold on
plot(h)
v=full(g.adjacency)/2;
%matrice periodica di periodo due perchè ci sono due primi vicini
sym(v/2)
sym((v)^1000)
%se elevo a una potenza molto alta tutti i contributi dovuti agli autovalori
%di modulo minore di 1 sono spariti, rimanendo solo quelli di modulo 1
% se voglio avere un equilibrio statistico devo sommare a w^1001/2
sym((v^1000+v^1001)/2)
%quello che succede è che la matrice diventa un proiettore sullo stato di
%equilibrio che è quello in cui la distribuzione è perfettamente uniforme
%(spiegato nelle note).
L1=full(g.laplacian)
%le considioni sono quelle libere: c'è differenza tra prima e ultima
%colonna e il balk
L2=full(h.laplacian)
%sia il laplaciano di g che il laplaciano di h sono matrici
%simmetriche, una con condizioni libere al bordo e l'altro con condizioni
%non libere;
%nel caso periodico le colonne del laplaciano sono simili , traslate di uno
%shift tra loro e con i punti agli estremi collegati.
%Per entrambi i laplaciani la somma su ogni colonna fa 0 (evidente dalla definizione)
%e questa è la ragione per cui c'è sempre l'autovettore costante come modulo
%0, cioè annichilito dal laplaciano: se applico il laplaciano a un vettore
%colonna costante questo fa 0, perchè devo fare la somma sulla riga in quel caso ma la
%matrice è simmetrica e quindi le somme su righe e colonne sono 0 per
%definizione. 
x=linspace(-5,5,256);
y=exp(-x.^2/2); fy=fft(y);
figure;
plot(x,y)
%c'è differenza a considerare la gaussiana tra -inf e +inf?
y(1)
%per calcolare i modi di fourier mi basta il primo dominio, ma posso
%prendere l'antitrasformata di fourier che mi da f(x) dati i modi di
%fourier, per ogni x anche oer x>5.
%attenzione: se ho una sola gaussiana, devo fare la trasformata di fourier e
%non la serie di fourier;
%per una funzione come la gaussiana molto ben localizzata, la trasformata
%di fourier e la serie di fourier tendono a essere la stessa cosa, sono la
%stessa cosa a meno di piccole correzioni.
figure;
plot(real(fftshift(fy)))
figure;
plot((-1).^(0:255).*real(fftshift(fy)),".-")
figure;
plot((-1).^(0:255).*imag(fftshift(fy)),".-")
%plottando la parte immaginaria cosa accade ?
%diverso da 0 e non è roundoff error, se lo fosse signica che fft fa schifo
%che non è vero. mi devo aspettare una parte immaginaria con un roundoff
%error di 10^-12 ecc.



