function [M,cp,kh2] = multidiagL(p)
% [M,cp,kh2] = multidiagL(p) returns the matrix M relevant for the 
% determination of the multidiagonal 1D Laplacian, the corresponding
% coefficients cp (in symbolic form) and the symfun kh2(k) which
% characterizes the spectrum.
