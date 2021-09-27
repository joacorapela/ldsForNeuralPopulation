[X,Y] = meshgrid(-100:20:100,-100:20:100);
 figure;
 %figure;quiver(X,Y,-1*X,-3*Y,1.5,'Color',[0.5,0.5,0.5])
 
%[a,b]=eig([-1.6,-1;-1,-1.6])
%hold on; line([-20,20],[-20,20],'Color','k','LineStyle','--','LineWidth',2)
hold on; quiver([-120],[-120],240,240,1.1,'Color',[0,0,0],'LineWidth',2)
%hold on; line([-20,20],[20,-20],'Color','k','LineStyle','--','LineWidth',2)
hold on; quiver([-120],[120],240,-240,1.1,'Color',[0.,0,0],'LineWidth',2)

%hold on;quiver(X,Y,-1.6*X-1*Y,-1*X-1.6*Y,1.6,'Color',[0.5,0.5,0.5],'LineWidth',0.5);
hold on;quiver(X,Y,-1.6*X-1*Y,-1*X-1.6*Y,1.2,'Color',[0.5,0.5,0.5],'LineWidth',1.3);
%%
[X,Y] = meshgrid(-100:20:100,-100:20:100);
figure;
hold on;quiver(X,Y,-1*X,-4*Y,0.8,'Color',[0.5,0.5,0.5],'LineWidth',1);
hold on; line([0,0],[-100,100],'Color','k','LineStyle','-','LineWidth',0.7)
hold on; line([-100,100],[0 0],'Color','k','LineStyle','-','LineWidth',0.7)

f=gcf;
f.Children.XTick = [];
f.Children.YTick = [];
set(f,'Color','w')
%%
figure;plot(0:0.01:1,exp(-5*(0:0.01:1)))
hold on;plot(0:0.01:1,exp(-20*(0:0.01:1)))
f=gcf;
f.Children.XTick = [];
f.Children.YTick = [];
set(f,'Color','w')
box('off')