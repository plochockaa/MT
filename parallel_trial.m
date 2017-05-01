clear all; close all; clc;

%> INPUTS : Nactin, cell_ecc
Nactin = 130;  gamma = [0 80 110 0]/180*pi;  ecc_cell = 0.7;  a = 60;
Nsegments = 15;
N = 132;

%%% setup using bio results
theta_actin = (25.2*pi/180)*(ecc_cell==0.7) + (19.8*pi/180)*(ecc_cell==0.75) + ...
         (14.7*pi/180)*(ecc_cell==0.8) + (12.8*pi/180)*(ecc_cell==0.85) + ...
         (7.7*pi/180)*(ecc_cell==0.9) + (8.7*pi/180)*(ecc_cell==0.95); % angle of parallel actin

% OUTPUTS: site(actin) to add to the end of the site matrix!!



%%%%%%%%% FIND GAP BETWEEN ACTIN FILAMENTS FOR ACTIN NUMBER %%%%%%%%%%%%%%%

ecc_standard = 0.98;
[b_ind] = lengths_for_ecc(ecc_standard,N,gamma,a,Nsegments);
L_ind = [a b_ind  a b_ind]; 
[xseed1,seedside1,L1,gamma1,sidenode1] = nodes_initiation(N,Nsegments,L_ind,gamma);

length1 = sqrt((sidenode1(1,1)-sidenode1(3,1))^2 + (sidenode1(1,2)-sidenode1(3,2))^2);
length2 = sqrt((sidenode1(2,1)-sidenode1(4,1))^2 + (sidenode1(2,2)-sidenode1(4,2))^2);
Lcell = max([length1,length2,L1]);

inc = Lcell/(Nactin+1);


%%%%%%%% NOW WORK WITH THE CELL OF A GIVEN ECC %%%%%%%%%%%%%%%%%%%%%%%%%%%%
[b] = lengths_for_ecc(ecc_cell,N,gamma,a,Nsegments);
L     = [a b  a b]; 
[xseed,seedside,L,gamma,sidenode] = nodes_initiation(N,Nsegments,L,gamma);
gamma1 = [0,cumsum(gamma(2:end))];

X = [xseed(:,1,1)-mean(xseed(:,1,1)),xseed(:,1,2)-mean(xseed(:,1,2))]; [u,s,v] = svd(X);
mid_cell = [mean(xseed(:,1,1)),mean(xseed(:,1,2))];

theta_cell = atan(v(2,1)/v(1,1));

theta_ACT = theta_cell-(pi/2-theta_actin);
theta_ACT = mod(theta_ACT,pi);

% describe the line at the main axis 
mid_x = mid_cell(1) + 2*Lcell*cos(theta_cell+pi);
mid_y = mid_cell(2) + 2*Lcell*sin(theta_cell+pi);

mid_r = 4*Lcell;

%%%%%%%%%%% FIND POSITIONS OF SIDENODES XING MAIN AXIS %%%%%%%%%%%%%%%%%%%%

x1 = mid_x;
y1 = mid_y;
theta1 = theta_cell;

    for kk = 1:4;

        x2 = sidenode(kk,1);
        y2 = sidenode(kk,2);
        theta2 = theta_ACT;
        [r1,r2] = intersectionrs01(x1,y1,theta1,x2,y2,theta2);
        R(kk) = (r1<0)*0 + (r1>=0).*r1;

    end

    min_R = min(R); % see figure 1 for the idea behind main axis calcs
    max_R = max(R);
    
    % actual number of actin filaments based on current cell ecc and
    % interval from ecc = 0.98
    actual_Nactin = floor((max_R - min_R)/inc);
    site     = zeros(N+actual_Nactin,Nsegments,4);  
    seed_actin = zeros(actual_Nactin,2);
    
    for jj = 1:actual_Nactin;
        
        seed_actin(jj,:) = [mid_x+(min_R+(jj)*inc)*cos(theta_cell),mid_y+(min_R+(jj)*inc)*sin(theta_cell)];
        
    end
    
    %%%%%%%%%%%%%%% FIND SEEDS ON SIDE EDGES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % move all the cell positions outside of the cell to find seeds easier!
    
    seed_actin2 = [seed_actin(:,1)+500*cos(theta_ACT+pi) seed_actin(:,2)+500*sin(theta_ACT+pi)];
    
    if theta_ACT<pi/2;
        
        for jj = 1:actual_Nactin;
        
        x1(jj) = seed_actin2(jj,1);
        y1(jj) = seed_actin2(jj,2);
        theta1 = theta_ACT;
   
           for kk = 2:3;
    
            x2 = sidenode(kk,1);
            y2 = sidenode(kk,2);
            theta2 = gamma1(kk);
            [r1,r2] = intersectionrs01(x1(jj),y1(jj),theta1,x2,y2,theta2);
            r(kk) = (r1<1.e-4)*600 + (r1>=1.e-4).*r1;
    
           end
        
        site(N+jj,1,2) = x1(jj) + min(r([2 3]))*cos(theta_ACT);
        site(N+jj,1,3) = y1(jj) + min(r([2 3]))*sin(theta_ACT);
        site(N+jj,1,4) = theta_ACT + pi;
    
        end
        
    elseif theta_ACT>pi/2;
        
          for jj = 1:actual_Nactin;
        
        x1(jj) = seed_actin2(jj,1);
        y1(jj) = seed_actin2(jj,2);
        theta1 = theta_ACT;
   
           for kk = 3:4;
    
            x2 = sidenode(kk,1);
            y2 = sidenode(kk,2);
            theta2 = gamma1(kk);
            [r1,r2] = intersectionrs01(x1(jj),y1(jj),theta1,x2,y2,theta2);
            r(kk) = (r1<1.e-4)*600 + (r1>=1.e-4).*r1;
    
           end
        
        site(N+jj,1,2) = x1(jj) + min(r([4 3]))*cos(theta_ACT);
        site(N+jj,1,3) = y1(jj) + min(r([4 3]))*sin(theta_ACT);
        site(N+jj,1,4) = theta_ACT + pi;
    
        end
        
        
    end
     
    xseed(N+1:N+actual_Nactin,1,1) = site(N+1:N+actual_Nactin,1,2);
    xseed(N+1:N+actual_Nactin,1,2) = site(N+1:N+actual_Nactin,1,3);
   
%%%%%%%%%%%% actin filaments (FIND RADIUS, we have pos and theta %%%%%%%%%%


  for jj = 1:actual_Nactin;
      
        x1a = site(N+jj,1,2);
        y1a = site(N+jj,1,3);
        theta1 = site(N+jj,1,4);
        
         for kk = 1:4;
    
            x2 = sidenode(kk,1);
            y2 = sidenode(kk,2);
            theta2 = gamma1(kk);
            [r1,r2] = intersectionrs01(x1a,y1a,theta1,x2,y2,theta2);
            r(kk) = (r1<1.e-4)*600 + (r1>=1.e-4)*r1;
    
         end
        
        if min(r)<600;
            
            site(N+jj,1,1) = min(r);
            
        else
            
            site(N+jj,1,1) = 0;
            
        end


  end

    
%%%%%%%%%%%
clf
figure(1)
hold all
plot([sidenode(:,1);  sidenode(1,1)], [sidenode(:,2);  sidenode(1,2)],'Color',[1 0 1],'LineWidth',3); % plot the cell sides 
scatter(mid_cell(1),mid_cell(2),'gx');
line([mid_x mid_x+mid_r*cos(theta_cell)],[mid_y mid_y+mid_r*sin(theta_cell)]);
for jj = 1:4;
scatter(mid_x+R(jj)*cos(theta_cell),mid_y+R(jj)*sin(theta_cell),'bo');
end
scatter(seed_actin(:,1),seed_actin(:,2),'go');

figure(2);
hold all;
plot([sidenode(:,1);  sidenode(1,1)], [sidenode(:,2);  sidenode(1,2)],'Color',[1 0 1],'LineWidth',3); % plot the cell sides 
scatter(seed_actin2(:,1),seed_actin2(:,2),'go');
for kk = 1:length(seed_actin2(:,1));
line([seed_actin2(kk,1),seed_actin2(kk,1)+ 800*cos(theta_ACT)],[seed_actin2(kk,2),seed_actin2(kk,2) + 800*sin(theta_ACT)]);
end

figure(3);
hold all;
plot([sidenode(:,1);  sidenode(1,1)], [sidenode(:,2);  sidenode(1,2)],'Color',[1 0 1],'LineWidth',3); % plot the cell sides 
scatter(site(N+1:N+actual_Nactin,1,2),site(N+1:N+actual_Nactin,1,3),'gx');

% hold all
% plot([sidenode(:,1);  sidenode(1,1)], [sidenode(:,2);  sidenode(1,2)],'Color',[1 0 1],'LineWidth',3); % plot the cell sides 
% scatter(x1,y1,'bx');
% scatter(seed_x,seed_y,'go');
% %scatter(actin_x,actin_y,'rx');
% axis('tight');
% 
% figure(2)
% hold all
% plot([sidenode(:,1);  sidenode(1,1)], [sidenode(:,2);  sidenode(1,2)],'Color',[1 0 1],'LineWidth',3); % plot the cell sides 
% plot([seed_x, seed_x+cos(theta_ACT-pi)*actin_r']',[seed_y, seed_y+sin(theta_ACT-pi)*actin_r']','b-','LineWidth',2);
% 
% figure(3)
% hold all
% plot([sidenode(:,1);  sidenode(1,1)], [sidenode(:,2);  sidenode(1,2)],'Color',[1 0 1],'LineWidth',3); % plot the cell sides 
% plot([seed_x, seed_x+cos(theta_ACT-pi)*actin_r']',[seed_y, seed_y+sin(theta_ACT-pi)*actin_r']','b-','LineWidth',2);
% %axis([0 sidenode(3,1) 0 sidenode(3,2)]);
% 
figure(4)
hold all;
plot([sidenode(:,1);  sidenode(1,1)], [sidenode(:,2);  sidenode(1,2)],'Color',[1 0 1],'LineWidth',3); % plot the cell sides 
for a=N+1:actual_Nactin+N;
    xactin=site(a,1,1)*cos(site(a,1,4));
    yactin=site(a,1,1)*sin(site(a,1,4)); 
          hold on
    plot([0 xactin]+xseed(a,1,1),[0 yactin]+xseed(a,1,2),'Color','b','LineWidth',1.5);
    scatter(site(a,1,2),site(a,1,3),'rx');

end