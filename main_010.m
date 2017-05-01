%==========================================================================
% The main code for the MT dynamics. 
%==========================================================================

clear all; close all; clc;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng('shuffle');
N         =  132;     % # of MT  
Nsegments =   20;     % 25;  

dt = .01/2;           % plotting time step 
draw0screen1 = 0 ;    % if draw0screen1 = 0, then the code draws, otherwise the output is to screen

t_end       = 10;     % 
drawfinal   = 0; 
drawall     = 1; 

thetacrmt   = 0 *pi/180;  % >>> the angle for zipping -- if set (thetacrmt = 0), then will have straight MTs (the code runs fater then, can reduce the number of segments to 1 or 2)
thetacrwall = thetacrmt;   % >>> the angle for zipping along the wall
%ecc_array =  [0.7 0.75 0.8 0.85 0.9 0.95 0.98];
ecc_array = 0.95; 

half_width = 1.5;   

xseed1 = zeros(N,2);   % "seeds"= minus ends of the MTs

tic % start counting how much time the code takes 

%%%%%%%%%%%%%%%%%%%%%%% TYPE OF SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nactin = [0];

% Type of actin: 1 - random, 2 - parallel
type_act = 2; 

% if actin_collision = true >> this is the new actin rules (pcat and cross)
actin_collision = true; 

% if zipbulk = true then COLLAPSE/ZIP < theta_crit & if zipbulk = false then CROSS/ZIP/COLLAPSE
zipbulk = false; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dir_name = ['/Users/aleksandraplochocka/Desktop/CELL/'];

haha_array = []; % collect the output of the main code 

     
 for kk =1 :length(ecc_array), 


    ecc = ecc_array(kk) ;
    
    Pcat = 0.01;
    
    %______( Markov chain transition rates )_____________

    alphapr  = 4; 
    betapr   = 1;
    alpha    = 2*500;  
    beta     = 3.5*alpha; 
    avglengthmttheory =  beta.*(alpha + betapr) ./ (beta*betapr - alpha * alphapr); 
   
    casestudy = 1; 
    
    %_______________( initiating the cell )____________________________
    a   = 60; 
    
    gamma = [0   80 110 0  ]/180*pi;
    
    [b] = lengths_for_ecc(ecc,N,gamma,a,Nsegments);
    % >> these 'b' calcs work for ecc 0.7 and higher only!! << FIX THIS
    
    Lsize = max(10, 1/2* beta*(alpha + betapr)./(beta*betapr  - alpha*alphapr)); 
    L     = [a b  a b]; 
   
    
    timetodraw = 0; 
     
    Nsides   = length(L); 
    xseed    = zeros(N,2);   % N = number of MTs, xseed = seed positions of the N MTs
    seedsdie = []; 
    sidenode = []; 

    if Nactin>0;
        
        if type_act == 1;
            
            [xseed,seedside,L,gamma,sidenode,site] = nodes_initiation_rand(N,Nactin,Nsegments,L,gamma);
            
        elseif type_act == 2;
            
            % >> fix these so that the ones out of the boundary are not shown
          
            [xseed,seedside,L,gamma,sidenode,site] = nodes_initiation_parallel(N,Nactin,Nsegments,L,gamma,theta_actin);
        end
        
    else
        
    [xseed,seedside,L,gamma,sidenode] = nodes_initiation(N,Nsegments,L,gamma);
    site     = zeros(N+Nactin,Nsegments,4);         % initialize all the sites have length zero  
    
    end
    
    % finding the main axis of the cell; xseed1 array is introduced ONLY for this chunk of the code. 
    for ii = 1:N,  xseed1(ii,1) = xseed(ii,1,1);  xseed1(ii,2) = xseed(ii,1,2); end
    cc = mean(xseed1); 
    [U,S,V] = svd([xseed1(:,1) - cc(1), xseed1(:,2) - cc(2)],0);
    avguu                   = zeros(1000,1); 
    uustd                   = zeros(1000,1); 
    uustd0                  = zeros(1000,1); 
    drawingtimearray        = zeros(1000,1); 
    
    
    
    %>>>> write on the screen we start running now. <<<< 
    disp(['>>> NEW RUN <<< , ecc = ',num2str(ecc), ', L = ', num2str(L)]);

    
 
    for i = 1:1,  % set the number of runs here for statistics, e.g. [for i = 1:200,] 
       
         haha = proga_zippering2...
               (timetodraw, dt, L, gamma, Nsides, xseed, seedside, sidenode, ...
                N,Nsegments,alpha, alphapr, beta, betapr, t_end, drawfinal,...
                drawall, casestudy,thetacrwall,thetacrmt,ecc,avglengthmttheory,xseed1,Pcat,draw0screen1,half_width,dir_name);
         
         
       
    end
     
  'done'
  haha_array = [haha_array,haha];
  
 end


toc


