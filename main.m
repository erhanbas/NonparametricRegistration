% Rigid registration with mutual information using KDE
% Algorithm first estimates the kernel density, then calculates the MI,
% maximizes MI using direct search method
% Erhan Bas
% 02/25/10

clear all;close all;clc

% Rigid transformation model
% R(y-t) = x, where y is unregistered and x is registered image. R is
% defined counterclockwise
% For images rotatioin is defined wrto image center

% Load Data
% 1) sythetic rectangle
% 2) brain MRI
% 3) manually selected registration curves

dataType = 1;
if dataType == 1
    N1 = 200;
    x = [5 0;0 3]*rand(2,N1);%-2.5;
    theta = pi/9;
    tx = 8; ty = 4;
    R = [cos(theta) sin(theta);-sin(theta) cos(theta)];
    y = R'*x(:,1:1:end) + repmat([tx;ty],1,length(x));
    N2 = length(y);
    anisotropic = 1;
    wgt = [theta;tx;ty];
    
    figure,
    plot(x(1,:),x(2,:),'.',y(1,:),y(2,:),'r.'), legend('x','y')
  
elseif dataType ==2
    unregistered = imread('unreg.bmp');
    registeredOriginal = imread('reg.bmp');
    figure, subplot(211),imshow(unregistered),title('y = unreg'),subplot(212), imshow(registeredOriginal),title('x = reg');
    I1 = edge(rgb2gray(unregistered),'canny',.4);
    I2 = edge(rgb2gray(registeredOriginal),'canny',.4);
    figure, subplot(211),imshow(I1),title('y = unreg'),subplot(212), imshow(I2),title('x = reg');
    [y1 y2]=find(I1); y = [y2 y1]';
    [x1 x2]=find(I2); x = [x2 x1]';
    
    y = y(:,1:6:end);
    x = x(:,1:6:end);
    
    N2 = length(y);
    anisotropic = 1;
    
    thetagt = pi/9;
    Rgt = [cos(thetagt) sin(thetagt);-sin(thetagt) cos(thetagt)];
    tgt =(eye(2)-Rgt')*[128;128];
    tgtx = tgt(1);
    tgty = tgt(2);
    wgt = [thetagt;tgt];
    xgt = Rgt*(y-repmat(tgt,1,length(y)));
    
    figure,
    plot(x(1,:),-x(2,:),'.',y(1,:),-y(2,:),'r.'), legend('x','y')
    
    figure,
    plot(x(1,:),-x(2,:),'.',xgt(1,:),-xgt(2,:),'r.'), legend('x','xgt')
    
elseif dataType ==3
    KV = load('EdgeSelection/edgesKV');
    DRR = load('EdgeSelection/edgesDRR2');
    x = KV.edgesKV;
    y = DRR.y;
    N2 = length(y);
    anisotropic = 1;
    
    figure,
    plot(x(1,:),-x(2,:),'.',y(1,:),-y(2,:),'r.'), legend('x','y'),
    
    %     estimating initial values of R and t
    %  x =R(y-t)
    % 1) Least mean square solution : Ax = b => x = inv(A'A)A'b
    % [x1;x2][cos sin;-sin cos][y1;y2]-[cos sin;-sin cos][t1;t2]
    % A = [y(1,1) y(1,2) 1 0;
    %      y(1,2) -y(1,1) 0 1;
    %      y(2,1) y(2,2) 1 0;
    %      y(2,2) -y(2,1) 0 1]
    % b = [cos;sin;k(1);k(2)] where k = -Rt
    
    % 2) Assume rotation is 0, match center locations
    
    mx_y = mean(y,2)-mean(x,2);
    x_ = x-repmat(mean(x,2),1,length(x));
    y_ = y-repmat(mean(y,2),1,length(y));

    init3 = [0*(1+.1*randn()) mx_y(1)*(1+.001*randn()) mx_y(2)*(1+.001*randn())];
    
    figure,
    plot(x_(1,:),-x_(2,:),'.',y_(1,:),-y_(2,:),'r.'), legend('x','y'),title('mean subtracted')
end
figure,
plot(x(1,:),-x(2,:),'.',y(1,:),-y(2,:),'r.'), legend('x','y')
% matlabpool(4)
%%
% Algorithm part
close all
if (1)
    [sigma_optx,S_x] = fitkdeFast(x,anisotropic);
    [sigma_opty,S_y] = fitkdeFast(y,anisotropic);
    scale=4;
    %%
    tic
    myfun = @(v) kdeopt2(v,y,x,S_x,S_y,anisotropic,scale);
    switch dataType
        case 1
            init = [wgt(1)*(1+.1*randn()) wgt(2)*(1+.1*randn()) wgt(3)*(1+.1*randn())];
        case 2
            init = [wgt(1)*(1+.05*randn()) wgt(2)*(1+.1*randn()) wgt(3)*(1+.1*randn())];
        case 3
            init = init3;
    end
    w2 = fminsearch(myfun,init);
    toc
    
    figure,
    plot(x(1,:),x(2,:),'.',y(1,:),y(2,:),'r.'), legend('x','y'), title('original')
    w=w2;
    
    the = w(1);
    ttx = w(2);
    tty = w(3);
    Rest = [cos(the) sin(the);-sin(the) cos(the)];
    xest = (Rest*(y - repmat([ttx;tty],1,length(y))));
    
    figure,plot(xest(1,:),-xest(2,:),'r.',x(1,:),-x(2,:),'o'), legend('xest','x'), title('forward estimation')
    yest = (Rest'*x + repmat([ttx;tty],1,length(x)));
    figure,plot(yest(1,:),-yest(2,:),'.b',y(1,:),-y(2,:),'r.'), legend('yest','y'), title('backward estimation')
end




