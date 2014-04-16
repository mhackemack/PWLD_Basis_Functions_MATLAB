close all; clear all; clc;

% number of vertices
nv=8;

% generate random points
x=randn(nv,1);
y=randn(nv,1);

% re-arrange so that there are in ccw order
xmean=mean(x);
ymean=mean(y);
[dummy, i] = sort(atan2(y-ymean,x-xmean));
x = x(i);
y = y(i);

% plot it for the sake of it
plot(x([1:end 1]),y([1:end 1]),'+r-')
hold on
for i=1:length(x)
    str=sprintf('%d',i);
    text(x(i)*1.05,y(i)*1.05,str);
end

% generate pwld matrices
[M,K,f,g]=build_pwld_local_matrices(1:nv,[x y]);
[x y]
M
K
f
for i=1:nv
    for j=1:nv
        gx(i,j)=g(1,i,j);
        gy(i,j)=g(2,i,j);
    end
end
gx
gy
