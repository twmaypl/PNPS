
N=50; xmin=-1; xmax=1;
CGL_x = -cos(pi*(0:N)'/N); % Chebyshev-Gauss-Labatto Points
[LGL_x,w,P]=lglnodes(N);  % Legendre-Gauss_Lobatto
LGL_x=flipud(LGL_x);  % LGL_x: Legendre-Gauss-Lobatta grid points
w=flipud(w);          %     w: LGL quadrature weights
M = diag(w); Minv = inv(M);  % M: Mass matrix; Minv: M^{-1}

[dpy, ~] = lepoly(N,LGL_x);    % First derivative of legendre polynomials dy
[dpx, ~] = lepoly(N,CGL_x);    % First derivative of Chebyshev polynomials dx
lc_p = ((1+LGL_x) .* dpx(:))/(2*dpy(N+1)); % Legendre intepolation polynomial using Chebyshev points
lc_m = ((1-LGL_x) .* dpx(:))/(2*dpy(N+1)); % Legendre intepolation polynomial using Chebyshev points

h = LegIntPol(CGL_x,LGL_x,N); % Legendre Interpolation Polynomial
S = h';          % S is the lgl to cgl transform matrix 
T = inv(h');     % T is the cgl to lgl transform matrix
                      
% set up differentiation grid   points 
% DiffMat=collocD(LGL_x);         % D: Differentiation matrix
[DiffMat_cgl,~] = cheb(N); 
DiffMat_cgl = -DiffMat_cgl;% CGL Differentiation matrix

% setup grid points in physical space and compute transformation metrics
x = xmin + (xmax-xmin)/(CGL_x(N+1)-CGL_x(1)) * (CGL_x(1:N+1)+1);
% y = xmin + (xmax-xmin)/(LGL_x(N+1)-LGL_x(1)) * (LGL_x(1:N+1)+1);
dx_dxi = chebfft1(x,1); dxi_dx = 1./dx_dxi; Jacobian= dx_dxi;


