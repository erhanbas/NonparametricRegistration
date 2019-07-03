nsig=.3;
N1=200;t=linspace(-pi,pi,N1);x=sin(t);y=cos(t);
cluster1=[x;y]+nsig*rand(2,N1);
N1=200;t=linspace(-pi,pi,N1);x=sin(t)+1;y=-cos(t);
cluster2=[x;y]+nsig*rand(2,N1);
data=[cluster1 cluster2];
%%
nsig = 3;
N1=200;cluster1=nsig*rand(2,N1);
N1=200;tet=pi/6;cluster2=[cos(tet) -sin(tet);sin(tet) cos(tet)]*cluster1+3;

data=[cluster1 cluster2];

[dim N]=size(data);
figure,plot(data(1,:),data(2,:),'.'),axis equal;

% [K Kx_xi]=kernel_matrix(data1,data2,n,N1,N2,params)