clear all; close all; clc;

N = 132;
ecc = 0.95;
a   = 1; 
b   = a./sqrt(1 - ecc^2); 


Nactin = 100;
Nsegments = 12;

[b] = lengths_for_ecc(ecc,N,gamma,a);

L     = [a b  a b]*60; 

xseed  = zeros(N+Nactin,2); 

[xseed,seedside,L,gamma,sidenode] = nodes_initiation(N,Nsegments,Nactin,L,gamma);

gamma1 = [0,cumsum(gamma(2:end))];

xv = [sidenode(1,1);sidenode(2,1);sidenode(3,1);sidenode(4,1);sidenode(1,1)];
yv = [sidenode(1,2);sidenode(2,2);sidenode(3,2);sidenode(4,2);sidenode(1,2)];

min_x = max(xv);
max_x = max(xv);
min_y = max(yv);
max_y = max(yv);

pos_actin = zeros(Nactin,2);

while sum(pos_actin(:,1)==0)>0;
    
    % generate random points in rectangle with boundaries defined by outer
    % points of trapezium
    
    if min_x < 0;
        L = max_x + min_x;
        pos_x = rand*L + min_x;
    else
        pos_x = rand*max_x;
    end

    pos_y = rand*max_y;
    
    [in,on] = inpolygon(pos_x,pos_y,xv,yv);
    
    if in == 1;
        
        index = find(pos_actin(:,1)==0);
        pos_actin(index(1),:) = [pos_x pos_y];
        
    end
    
end

for jj = 1:Nactin;

        r = zeros(1,4);
        x1(jj) = pos_actin(jj,1);
        y1(jj) = pos_actin(jj,2);
        theta1(jj) = rand*2*pi;
        
        for kk = 1:4
            x2 = sidenode(kk,1);
            y2 = sidenode(kk,2);
            theta2 = gamma1(kk);
            [r1,r2] = intersectionrs01(x1(jj),y1(jj),theta1(jj),x2,y2,theta2);
            r(kk) = (r1<0)*600 + (r1>=0).*r1;
        end
        
       radius(jj) = min(r);
       x(jj) = radius(jj)*cos(theta1(jj));
       y(jj) = radius(jj)*sin(theta1(jj));
       


end

seeds_actin = [[pos_actin(:,1)+x'], [pos_actin(:,2)+y']];

for jj = 1:Nactin;
    

        r = zeros(1,4);
        x1(jj) = seeds_actin(jj,1);
        y1(jj) = seeds_actin(jj,2);
        theta(jj) = mod(theta1(jj)+pi,2*pi);
        
        for kk = 1:4
            x2 = sidenode(kk,1);
            y2 = sidenode(kk,2);
            theta2 = gamma1(kk);
            [r1,r2] = intersectionrs01(x1(jj),y1(jj),theta(jj),x2,y2,theta2);
            r_new(kk) = (r1<=1.e-5)*600 + (r1>1.e-5).*r1;
        end
        
       radius_new(jj) = min(r_new);
       x_new(jj) = radius_new(jj)*cos(theta(jj));
       y_new(jj) = radius_new(jj)*sin(theta(jj));
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    
end

xseed(1:Nactin,1,1) = seeds_actin(:,1);
xseed(1:Nactin,1,2) = seeds_actin(:,2);
site(1:Nactin,1,1) = theta;
site(1:Nactin,1,2) = seeds_actin(:,1);
site(1:Nactin,1,3) = seeds_actin(:,2);
site(1:Nactin,1,4) = radius_new;

%%%%%%%%%%%
clf
figure(1)
hold all
plot([sidenode(:,1);  sidenode(1,1)], [sidenode(:,2);  sidenode(1,2)],'Color',[1 0 1],'LineWidth',3); % plot the cell sides 
scatter(pos_actin(:,1),pos_actin(:,2),'bx','LineWidth',2);
plot([pos_actin(:,1) , pos_actin(:,1) + x']',...
[pos_actin(:,2) , pos_actin(:,2) + y']','g-','LineWidth',2);
scatter(seeds_actin(:,1),seeds_actin(:,2),'bo','LineWidth',2);
plot([seeds_actin(:,1) , seeds_actin(:,1) + x_new']',...
[seeds_actin(:,2) , seeds_actin(:,2) + y_new']','g--','LineWidth',2);


figure(2)
hold all
plot([sidenode(:,1);  sidenode(1,1)], [sidenode(:,2);  sidenode(1,2)],'Color',[1 0 1],'LineWidth',3); % plot the cell sides 
plot([seeds_actin(:,1) , seeds_actin(:,1) + x_new']',[seeds_actin(:,2) , seeds_actin(:,2) + y_new']','b-','LineWidth',2);
