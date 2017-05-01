clear all; close all; clc

N = 100; Nsegments = 2; a = 60; gamma = [0   80 110 0  ]/180*pi;

ecc_array = [0.98];




for jj = 1:length(ecc_array);
        ECC2 = ecc_array(jj)^2;

        change = zeros(1,1000000);

        t = 1;
        change(1) = 0;
        inc = 0.5;
        b   = a./sqrt(1 - ECC2); 
        B_og(jj) = b;
        L     = [a b  a b]; 

 
        [xseed,seedside,L,gamma,sidenode] = nodes_initiation(N,Nsegments,L,gamma);

        X = [xseed(:,1,1)-mean(xseed(:,1,1)), xseed(:,1,2)-mean(xseed(:,1,2))];
        [u,s,v] = svd(X);
        
        
        ecc_new2(t) =  1 - (s(2,2)/s(1,1))^2;


        while (abs(ecc_new2(t)-ECC2)>eps2) && (t < t_max);
             
            
            if ecc_new2(t) > ECC2;
                t = t+1;
                change(t) = -1;
                b = b - inc;

            else
                t = t+1;
                change(t) = 1;
                b = b + inc;
            end 
            

            L     = [a b a b]; 

            [xseed,seedside,L,gamma,sidenode] = nodes_initiation(N,Nsegments,L,gamma);

            X = [xseed(:,1,1), xseed(:,1,2)];
            X = [xseed(:,1,1)-mean(xseed(:,1,1)), xseed(:,1,2)-mean(xseed(:,1,2))];
            [u,s,v] = svd(X);
            ecc_new2(t) =  1- (s(2,2)/s(1,1))^2;
            
            figure(1); 
                subplot(221); 
                    plot(1:t,sqrt(ecc_new2(1:t))); 
                subplot(222);
                    plot(t,s(1,1),'o','markersize',10); hold on; 
                    plot(t,s(2,2),'d','markersize',10); hold on; 
                  
                subplot(223); 
                    line([0,v(1,1)],[0,v(1,2)],'Color','blue'); hold on; 
                    line([0,v(2,1)],[0,v(2,2)],'Color','magenta'); hold on; 
                    axis([-2 2 -2 2])
                    axis equal
                    grid on
                    
                subplot(224);
                    plot(xseed(:,1,1),xseed(:,1,2),'o','markersize',2);
%                      axis([-80 80 -80 80])
%                     pause(.1);

            if change(t)~=change(t-1);
                inc = inc/2;
            end

 

        B(jj) = b;

        
        end
        
        
end
        
% for kk = 1:length(ecc_array);
%     
%     L     = [a B(kk)  a B(kk)]; 
%      
% [xseed,seedside,L,gamma,sidenode] = nodes_initiation(N,Nsegments,L,gamma);
% 
% X = [xseed(:,1,1)-mean(xseed(:,1,1)),xseed(:,1,2)-mean(xseed(:,1,2))]; [u,s,v] = svd(X);
% theta_cell(kk) = atan(v(2,1)/v(1,1));
% 
% end

% theta_cell = theta_cell*180/pi;
% 
% dir_name = ['/Users/aleksandraplochocka/Desktop/MTSD/'];
% 
% 
% filename_mtsd = [dir_name,'cell_theta.csv'];
% csvwrite(filename_mtsd, theta_cell);

