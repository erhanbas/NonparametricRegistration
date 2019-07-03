function [sigma_opt,cov_x] = fitkde(x,anisotropic)
% Uses Fibonacci search to identify optimal kernel size for variable-width KDE
figure(32),clf,
[D,L] = size(x);        % determine data dimensionality and sample size
if anisotropic == 1
    % Use local k-NN covariances around x_i
    k = max(ceil(sqrt(L)),D+1); PairDist = squareform(pdist(x'));
    for i = 1:L
        [dummy,indknn] = sort(PairDist(i,:),'descend'); xiknn = x(:,indknn(1:k));
        temp = repmat(x(:,i),1,k)-xiknn; cov_x{i} = temp*temp'/k;
    end
elseif anisotropic == 0
    cov_x = 0;
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
        Hs(length(sigmaold)+s) = EntropyEst(x,cov_x,sigma(length(sigmaold)+s),anisotropic);
    end,
    figure(32), plot(sigma(length(sigmaold)+1:end),Hs(length(sigmaold)+1:end),'.','color',rand(3,1)),hold on,drawnow,
    [sigma,ind] = sort(sigma,'ascend'); Hs = Hs(ind); [dummy,idx] = min(Hs);
    if (idx==1)|(idx==length(sigma)),
        flag = 0, disp('minimum outside initial interval');
    else,
        startpoint = log10(sigma(idx-1)); endpoint = log10(sigma(idx+1));
    end,
    if count > countmax, flag = 0; end,
    count = count + 1; sigmaold = sigma;
end
[dummy,tmpind] = min(Hs); disp('optimum sigma found...'), sigma_opt = sigma(tmpind),

%--------------------------------------------
function H = EntropyEst(x,cov_x,sigma,anisotropic)
[n,L] = size(x); sigma2 = sigma^2;
if anisotropic == 1
    S = zeros(n*L,n*L); C = zeros(1,L);
    for i = 1:L
        S((i-1)*n+1:i*n,(i-1)*n+1:i*n) = inv(sigma2*cov_x{i});
        C(i) = 1/((2*pi)^(n/2)*det(sigma2*cov_x{i})^(1/2));
    end
    for j = 1:L
        xjmx = repmat(x(:,j),1,L)-x; indj = [1:L]; indj(j) = [];
        tempj = reshape(S*reshape(xjmx,n*L,1),n,L);
        p(j) = sum(C(indj).*exp(-0.5*sum(xjmx(:,indj).*tempj(:,indj),1)))/(L-1);
    end
    H = -mean(log(p));
elseif anisotropic == 0
    C = 1/(2*pi*sigma2)^(n/2);
    for j = 1:L
        xjmx = repmat(x(:,j),1,L)-x; xjmx(:,j) = [];
        p(j) = C*sum(exp(-0.5*sum(xjmx.^2,1)/sigma2))/(L-1);
    end
    H = -mean(log(p));
end
%--------------------------------------------