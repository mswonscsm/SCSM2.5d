%Seismogram data is recorded in rec_real_xx.out
%rec_real_xx is matrix of [14,N-Timestep]
%1st row: time
%2nd row: water velocity x
%3rd row: water velocity y
%4th row: water velocity z
%5th row: water pressure
%6th row: solid velocity x
%7th row: solid velocity y
%8th row: solid velocity z
%9th row: solid stress xx
%10th row: solid stress xy
%11th row: solid stress yy
%12th row: solid stress xz
%13th row: solid stress yz
%14th row: solid stress zz

figure
t=tiledlayout("vertical");
t.TileSpacing = 'none';
numfiles = 7;
set(gcf,'Color','w');

for k = 2:numfiles
  myfilename = sprintf('rec_real_%d.out', k);
  rec = importdata(myfilename);

  nexttile(t)
  hold on
  plot(rec(1,:),rec(5,:) ,'k','LineWidth',1.5)
  
  set(gca,'YTick',[])
  if(k<numfiles)
  set(gca,'XTick',[])
  set(gca,'XColor','none')
  end
  pbaspect([5 1 1])
  
  if(k==2)
  %legend boxoff
  set(gca,'Fontsize',14)
  end
end
set(gca,'Fontsize',18)



