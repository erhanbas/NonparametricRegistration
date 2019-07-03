function p = kdeopt(v,y,x,Sx,Sy,anisotropic)
theta = v(1);
t = v(2);
R = [cos(theta) sin(theta);-sin(theta) cos(theta)];
[ny,Ly] = size(y); [nx,Lx] = size(x); 
if ny~=nx, error('Dimension mismatch!'), end, n = nx;
if anisotropic == 1
    DS = zeros(Lx,Ly);
    S = zeros(Lx,Ly,4);
    for i = 1:Lx
        Sx_=Sx((i-1)*nx+1:i*nx,(i-1)*nx+1:i*nx);
        for j = 1:Ly
            Sy_ = Sy((j-1)*ny+1:j*ny,(j-1)*ny+1:j*ny);
            S_ = (R'*Sx_*R + Sy_);
            DS(i,j) = 1/((2*pi)^(ny/2)*det(S_)^(1/2));
            S_ = inv(S_);
            S(i,j,:) = S_(:);
        end
    end
    K = zeros(Ly,1);
    for j = 1:Ly
        %         if mod(j,100)==0, [j,Ly], end,
        
        yjmx = (R'*x - t) -repmat(y(:,j),1,Lx) ;
        Yjmx = [yjmx(1,:).^2;yjmx(1,:).*yjmx(2,:);yjmx(1,:).*yjmx(2,:);yjmx(2,:).^2];
        K(j) = mean(DS(:,j).*exp(-.5*sum(squeeze(S(:,j,:))'.*Yjmx,1))');
    end
    p = -sum(K);
elseif anisotropic == 0
    C = 1/(2*pi*sigma2)^(n/2);
    for j = 1:Ly
        yjmx = repmat(y(:,j),1,Lx)-x;
        p(j) = C*mean(exp(-0.5*sum(yjmx.^2,1)/sigma2));
    end
end

