%Solid velocity x : fort.101, 102, 103...
%Solid velocity z : fort.301, 302, 303...
%Water pressure : fort.401, 402, 403...

load Xgrid.out
load Zgrid.out

xmin=min(Xgrid);xmax=max(Xgrid);
zmin=min(min(Zgrid));
zmax=max(max(Zgrid));

%%%%%%%%%%%%%
check=load ('fort.402');

%%%%%%%%%%%%%

nx=size(check,1);
nz=size(check,2);
mx=max(max(check));
mn=min(min(check));

for i=1:nx
    for j=1:nz
        if check(i,j)<0
            check(i,j)=-check(i,j)/mn;
        elseif check(i,j)>0
            check(i,j)=check(i,j)/mx;
        end
     end
end

%%
x1=zeros(size(Zgrid));
for i=1:nz
    x1(:,i)=Xgrid;
end

figure(1),surface(x1,Zgrid,check,'EdgeColor','none')
colormap redbluecmap;
set(gcf,'Color', 'w');

set(gca,'XColor','none')
set(gca,'YColor','none')

daspect([1 1 1])
