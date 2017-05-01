% clear all; close all; clc;
% 
% N = 132;
% ecc = 0.95;
% a   = 1; 
% b   = a./sqrt(1 - ecc^2); 
% L     = [a b  a b]; 
% gamma = [0   80 110 0  ]/180*pi;
% Nactin = 150;
% Nsegments = 12;
% xseed  = zeros(N+Nactin,2); 
% gamma1 = [0,cumsum(gamma(2:end))];
% 
% [xseed,seedside,L,gamma,sidenode] = nodes_initiation(N,Nsegments,Nactin,L,gamma);

function [xseed, site] = nodes_random(xseed,seedside,L,gamma,sidenode,Nactin,N,Nsegments);
% This function outputs the site and xseed information for the randomized
% actin array.

 
 site     = zeros(N+Nactin,Nsegments,4);         % initialize all the sites have length zero  
 gamma1 = [0,cumsum(gamma(2:end))];
 
    for k=1:Nactin;
       
        %setup actin middle seed inside region of the cell
        if sidenode(3,2)<=sidenode(4,2);                                              
            midseed(k,2)=(sidenode(3,2)*0.5)*rand+(sidenode(3,2)/4);
        else
            midseed(k,2)=(sidenode(4,2)*0.5)*rand+(sidenode(4,2)/4);
        end
        
        midseed(k,1)= (cos(gamma(2))*midseed(k,2)) + rand*L(1)*0.5 + (L(1)/4);
  
    end
    
 %%%% Trial plot of this function    
%  plot([sidenode(:,1);  sidenode(1,1)], [sidenode(:,2);  sidenode(1,2)],'Color',[1 0 1],'LineWidth',3, 'LineSmoothing','on'); % plot the cell sides 
%  hold on; 

for a=N+1:Nactin+N;
    angle=rand*2*pi;
    r=zeros(1,4);
   
    %%%%% compare the middle actin seed to every side
        x1 = midseed(a-N,1); 
        y1 = midseed(a-N,2);
        theta1 = angle;
        
        for kk = 1:4
            x2 = sidenode(kk,1);
            y2 = sidenode(kk,2);
            theta2 = gamma1(kk);
            [r1,r2] = intersectionrs01(x1,y1,theta1,x2,y2,theta2);
            r(kk) = (r1<0)*600 + (r1>=0).*r1;
        end
    radius = min(r);
    
    xactin=radius*cos(angle+pi);  % from the initial angle add pi to rotate the position
    yactin=radius*sin(angle+pi);
    
    xseed(a,1,1) = (radius*cos(angle))+midseed(a-N,1); 
    xseed(a,1,2) = (radius*sin(angle))+midseed(a-N,2);
    
        x1 = xseed(a,1,1);   %check the maximum radius for actin, using position along cell boundary as seeds 
        y1 = xseed(a,1,2);
        theta1 = angle+pi; 
        m=zeros(1,4);
        

        for jj = 1:4        % update radius to for the seed, covering the length of the cell
            x2 = sidenode(jj,1);
            y2 = sidenode(jj,2);
            theta2 = gamma1(jj);
            [r1,r2] = intersectionrs01(x1,y1,theta1,x2,y2,theta2);
            m(jj) = (r1<0)*600 + (r1>=0).*r1;
        end
         radius = min(m(m>0.5));
         
    xactin=radius*cos(angle+pi);  
    yactin=radius*sin(angle+pi);
         
%     plot([0 xactin]+xseed(a,1,1),[0 yactin]+xseed(a,1,2),'Color','b','LineWidth',1.5, 'LineSmoothing','on');
%  hold on

site(a,1,1) = radius;
site(a,1,2) = xseed(a,1,1);
site(a,1,3) = xseed(a,1,2);
site(a,1,4) = theta1;

end

% hold off

end