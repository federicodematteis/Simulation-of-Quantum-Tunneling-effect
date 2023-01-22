function f = isfft(tf)
% inverse symmetric fast-fourier
% tf = isfft(f)
% acts columnwise unless f is a row vector

[N,nc] = size(tf);
isrow = 0;
if N == 1
    isrow = 1;
    tf = tf.';
    N = nc;
    nc = 1 ;
end

n = floor(N/2);
tf = tf([n+1:N,1:n],:);

n = floor((N-1)/2);
f = ifft(exp(-1i*pi*(1-1/N)*[0:n,n+1-N:-1]')*ones(1,nc) .* tf);

if isrow
    f = f.';
end
