clear all; clf;
N = 132; Nsegments = 1; Nactin = 55; 
ecc = 0.95;
a   = 1; 
b   = a./sqrt(1 - ecc^2); 


L     = [a b  a b]; 
gamma = [0   80 110 0  ]/180*pi;

for ii = 1:150
[xseed,seedside,L,gamma,sidenode] = nodes_initiation(N,Nsegments,Nactin,L,gamma);
[xseed, site] = nodes_random(xseed,seedside,L,gamma,sidenode,Nactin,N,Nsegments);


 figure(1);
 plot([sidenode(:,1);  sidenode(1,1)], [sidenode(:,2);  sidenode(1,2)],'Color',[1 0 1],'LineWidth',3, 'LineSmoothing','on'); % plot the cell sides 
 hold on; 

for a=N+1:Nactin+N;
    xactin=site(a,1,1)*cos(site(a,1,4));
    yactin=site(a,1,1)*sin(site(a,1,4)); 
          
    plot([0 xactin]+xseed(a,1,1),[0 yactin]+xseed(a,1,2),'Color','g','LineWidth',1, 'LineSmoothing','on');
    axis equal
    hold on
end
hold off
pause(1)
end 