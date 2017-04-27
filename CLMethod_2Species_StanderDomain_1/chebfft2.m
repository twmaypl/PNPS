function d2u = chebfft2(k,u,n)
% DiffMatt * K * DiffMatt * U 
% d/dx (k du/dx)

du = chebfft1(u,n);
d2u = chebfft1(k*du,n);
end