function p = kdeopt2(v,y,x,Sx,Sy,anisotropic,scale)
theta = v(1);
tx = v(2);
ty = v(3);
R = [cos(theta) sin(theta);-sin(theta) cos(theta)];
[ny,Ly] = size(y); [nx,Lx] = size(x); 
if ny~=nx, error('Dimension mismatch!'), end, n = nx;
if anisotropic == 1
    DS = zeros(Lx,Ly);
    S = zeros(Lx,Ly,nx*ny);
    parfor i = 1:Lx
        Sx_=Sx{i}*scale;
        DS_ = zeros(1,Ly);
        S2 = zeros(Ly,ny^2);
        for j = 1:Ly
            Sy_ = Sy{j}*scale;
            S_ = (R'*Sx_*R + Sy_);
            DS_(j) = 1/((2*pi)^(ny/2)*det(S_)^(1/2));
            S_ = inv(S_);
            S2(j,:) = S_(:);
        end
        DS(i,:) = DS_;
        S(i,:,:) = S2;
    end
    K = zeros(Ly,1);
    parfor j = 1:Ly
        %         if mod(j,100)==0, [j,Ly], end,
        yjmx = (R'*x + repmat([tx;ty],1,length(x))) -repmat(y(:,j),1,Lx) ;
        Yjmx = [yjmx(1,:).^2;yjmx(1,:).*yjmx(2,:);yjmx(1,:).*yjmx(2,:);yjmx(2,:).^2];
        K(j) = mean(DS(:,j).*exp(-.5*sum(squeeze(S(:,j,:))'.*Yjmx,1))');
    end
    p = -sum(K);
elseif anisotropic == 0
    C = 1/(2*pi*sigma2)^(n/2);
    p=zeros(1,Ly);
    for j = 1:Ly
        yjmx = repmat(y(:,j),1,Lx)-x;
        p(j) = C*mean(exp(-0.5*sum(yjmx.^2,1)/sigma2));
    end
end

