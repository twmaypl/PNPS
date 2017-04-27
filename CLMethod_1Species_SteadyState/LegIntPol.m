function h = LegIntPol(x,y,N)
% output the Legendre-Lagrange interpolation polynomial matrix at point x
% input y is a length N vector of legendre collocation points
% input x is a length N vector of other points
% l_i(x) = (1-X^2)P_N'(x) /(N(N+1)P_N(y_i)(x-y_i)

[~,P_y] = lepoly(N,y); % P_N'(y) (not output) & P_N(y)
[dP_x, ~] = lepoly(N,x); % P_N'(x) & P_(x) (not output)
h = zeros(N+1,N+1);
for k=1:N+1
    for j=1:N+1
        if(x(k) == y(j))
            h(j,k) = 1;
        else
            h(j,k) = -(1-x(k)^2) * dP_x(k) / (N*(N+1)*P_y(j)*(x(k)-y(j)));
%             if (h(j,k) < 1e-10)
%                 h(j,k) = 0;
%             end
        end
    end
end
