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

numfiles = 114;
load rec_real_1.out
nt=size(rec_real_1,2);
mydata = zeros(nt,numfiles-1);

x=[1:numfiles-1];
y=rec_real_1(1,:);

for k = 2:numfiles
  myfilename = sprintf('rec_real_%d.out', k);
  mydata0 = importdata(myfilename);
  mydata1=mydata0(5,:);
  mydata(:,k-1) = mydata1;
end

%%
%plotimage is the function of external toolbox. you can download below link 
%https://www.crewes.org/ResearchLinks/FreeSoftware/
plotimage(x,y(:),mydata(:,:))
xlabel('Line Trace')
ylabel('Time (s)')
colorbar

