function [sigma_opt,S_x] = fitkdeFast(x,anisotropic)
% Uses Fibonacci search to identify optimal kernel size for variable-width KDE
figure(32),clf,
[D,L] = size(x);        % determine data dimensionality and sample size
if anisotropic == 1
    % Use local k-NN covariances around x_i
    k = max(ceil(sqrt(L)),D+1); PairDist = squareform(pdist(x'));
    for i = 1:L
        [~,indknn] = sort(PairDist(i,:),'descend'); xiknn = x(:,indknn(1:k));
        temp = repmat(x(:,i),1,k)-xiknn; S_x{i} = inv(temp*temp'/k);
        C(i) = 1/((2*pi)^(D/2)*det(temp*temp'/k)^(1/2));
    end
elseif anisotropic == 0
    S_x = 0;
end
% Find the optimal kernel width using Fibonacci/Golden-Selection search in 10.^[-2,0.5]
F = [1 1 2 3]; count = length(F)+1; countmax = length(F)+2;
startpoint = -3; endpoint = 1; flag = 1; sigmaold = []; Hs = [];
while flag
    F(count) = F(count-1) + F(count-2); interval = (endpoint-startpoint)/F(count);
    linspc = [startpoint:interval:endpoint]; newind = [2:length(linspc)-1];
    sigma = [sigmaold,10.^linspc(newind)];
    for s = 1:length(newind),
        [s,length(newind)],
        Hs(length(sigmaold)+s) = EntropyEst(x,S_x,C,sigma(length(sigmaold)+s),anisotropic);
    end,
    figure(32), plot(sigma(length(sigmaold)+1:end),Hs(length(sigmaold)+1:end),'.','color',rand(3,1)),hold on drawnow,
    [sigma,ind] = sort(sigma,'ascend'); Hs = Hs(ind); [~,idx] = min(Hs);
    if (idx==1)|(idx==length(sigma)),
        flag = 0, disp('minimum outside initial interval');
    else,
        startpoint = log10(sigma(idx-1)); endpoint = log10(sigma(idx+1));
    end,
    if count > countmax, flag = 0; end,
    count = count + 1; sigmaold = sigma;
end
[~,tmpind] = min(Hs); disp('optimum sigma found...'), sigma_opt = sigma(tmpind),
for i = 1:L
    S_x{i} = S_x{i}/sigma_opt^2;
end

%--------------------------------------------
function H = EntropyEst(x,S_x,C,sigma,anisotropic)
[n,L] = size(x); sigma2 = sigma^2;
if anisotropic == 1
    S = zeros(n^2,L); %C = zeros(1,L);
    for i = 1:L
        temp = S_x{i}/sigma2;
        S(:,i) = temp(:);
        %         C(i) = 1/((2*pi)^(n/2)*det(sigma2*cov_x{i})^(1/2));
    end
    
    p = zeros(1,L);
%     tic
    for j = 1:L
        xjmx = repmat(x(:,j),1,L)-x; indj = [1:L]; indj(j) = [];
        Xjmx = [xjmx(1,:).^2;xjmx(1,:).*xjmx(2,:);xjmx(1,:).*xjmx(2,:);xjmx(2,:).^2];
        temp = (C/sigma2).*exp(-.5*sum(S.*Xjmx,1));
        p(j) = sum(temp(indj))/(L-1);
        %         tempj = reshape(S*reshape(xjmx,n*L,1),n,L);
        %         p(j) = sum((C(indj)/sigma2).*exp(-0.5*sum(xjmx(:,indj).*tempj(:,indj),1)))/(L-1);
    end
%     toc
     H = -mean(log(p));
   
%     tic
%     S = zeros(n*L,n*L); %C = zeros(1,L);
%     for i = 1:L
%         S((i-1)*n+1:i*n,(i-1)*n+1:i*n) = S_x{i}/sigma2;
%         %         C(i) = 1/((2*pi)^(n/2)*det(sigma2*cov_x{i})^(1/2));
%     end
%     toc
%     p2 = zeros(1,L);
%     tic
%     for j = 1:L
%         xjmx = repmat(x(:,j),1,L)-x; indj = [1:L]; indj(j) = [];
%         tempj = reshape(S*reshape(xjmx,n*L,1),n,L);
%         p2(j) = sum((C(indj)/sigma2).*exp(-0.5*sum(xjmx(:,indj).*tempj(:,indj),1)))/(L-1);
%     end
%     toc
%     H = -mean(log(p));
elseif anisotropic == 0
    C = 1/(2*pi*sigma2)^(n/2);
    for j = 1:L
        xjmx = repmat(x(:,j),1,L)-x; xjmx(:,j) = [];
        p(j) = C*sum(exp(-0.5*sum(xjmx.^2,1)/sigma2))/(L-1);
    end
    H = -mean(log(p));
end
%--------------------------------------------