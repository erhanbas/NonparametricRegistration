function p = kdeA(y,x,sigma,cov_x,anisotropic)
% Evaluates the KDE at y specified by the other inputs
[ny,Ly] = size(y); [nx,Lx] = size(x); sigma2 = sigma^2;
if ny~=nx, error('Dimension mismatch!'), end, n = nx;
if anisotropic == 1
    S = zeros(n*Lx,n*Lx); C = zeros(1,Lx);
    for i = 1:Lx
        S((i-1)*n+1:i*n,(i-1)*n+1:i*n) = inv(sigma2*cov_x{i});
        C(i) = 1/((2*pi)^(n/2)*det(sigma2*cov_x{i})^(1/2));
    end
    p = zeros(1,Ly);
    for j = 1:Ly
        if mod(j,100)==0, [j,Ly], end,
        yjmx = repmat(y(:,j),1,Lx)-x; tempj = reshape(S*reshape(yjmx,n*Lx,1),n,Lx);
        p(j) = mean(C.*exp(-0.5*sum(yjmx.*tempj,1)));
    end
elseif anisotropic == 0
    C = 1/(2*pi*sigma2)^(n/2);
    for j = 1:Ly
        yjmx = repmat(y(:,j),1,Lx)-x;
        p(j) = C*mean(exp(-0.5*sum(yjmx.^2,1)/sigma2));
    end
end

