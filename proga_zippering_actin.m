%==========================================================================
% Main zippering code (with lots of redundancies...) 
%
%
% notes: 
% - Minimum distance btw the seeds is 3
%
%   N         = number of active sites
%   alpha, alpha', beta, beta' = rates
%   t_end     = max length of the computationcl
%   drawfinal = (0/1) = draw the final version
%   drawall   = (0/1) = draw at every step.
%   casestudy = the number of the case scenario that MTs undertake; 
%               e.g. case = 1, then MT loses a cap upon bumping into the
%               wall.  
%   xseed1    = 
%   xseed     = 
%
%__________________________________________________________________________
%    site = [r, x, y, theta]; 
%            1  2   3  4       
%
% notes
%   r     - length of the MT
%   x     - x-unitary direction of the MT 
%   y     - y-unitary direction of the MT  
%   theta - angle wrt to the horizontal
%   x_s   - x of the seed (not along the MT, properly in 2D)
%   y_s   - y of the seed (not along the MT, properly in 2D)
%   side # -- number of the cell side that the MT belongs to 
%   
%
%==========================================================================
%
function qwerqwer  = proga_zippering_actin...
               (site, timetodraw, dt, L, gamma, Nsides, xseed, seedside, sidenode, ...
                N,Nsegments,Nactin,alpha, alphapr, beta, betapr, t_end, drawfinal,drawall, casestudy,...
                thetacrwall,thetacrmt,ecc,avglengthmttheory,xseed1,Pcat,draw0screen1,half_width,dir_name,testnumber,lines_measure);
            


cycles   = 0;                            % number of cycles that has started
state    = [alphapr, alpha,  alphapr;    % the first column is all the possible states when the cap is off
               beta, betapr, 0     ];    % the second column is all the possible states when the cap is on. 
                                         % the third column is the rate of exiting the zero length situation    
                                         

segm     = zeros(N+Nactin,1); 
segm(N+1:N+Nactin)=1;  %all Actin have one segment 

capyesno = zeros(N,1); 
bigstate = zeros((2*N),1); 
for i=1:N,  bigstate([2*i-1, 2*i]) = state(:,3); end % initally all the MT sites are in the B0 state, behaving like GDP, i.e. no cap.
time     = 0;
tindex   = 0; 

for ii = 1:N,  xseed1(ii,1) = xseed(ii,1,1);  xseed1(ii,2) = xseed(ii,1,2); end
cc = mean(xseed1); 
[U,S,V]          = svd([xseed1(:,1) - cc(1), xseed1(:,2) - cc(2)],0);
avguu            = zeros(1000,1); 
uustd            = zeros(1000,1); 
uustd0           = zeros(1000,1); 
drawingtimearray = zeros(1000,1); 
avglengthmt      = zeros(1000,1); 
mtdensity        = zeros(1000,1); 
fracnonzeroMTs   = zeros(1000,1); 
cellarea         = L(1)*L(2)*sin(pi - gamma(2))/2 + L(3)*L(4)*sin(pi - gamma(4))/2; 
            
theta_hist     = zeros(180,1);              % the histogram of the angles of all the segments 
theta_hist_avg = zeros(180,1);              % the averaged over time histogram of the angles of all the segments 



%_______________( initializing the intersection matrix)____________________
% Have two arrays: 
% - intsctysno = if the intersection happened between the i-th and j-th MT, 
% - intsctdist = how long to the intersection.

intsctinfty    = 600;                               % if there is no intersection, set the distance to intersection to be a big number  = 5 exp. length of  MT. 
intsctdist     = ones ((N+Nactin)*Nsegments,(N+Nactin)*Nsegments)*intsctinfty;  % matrix of intersections A(i,j) = if (i) intersects (j), distance for (i) to the intersection. 
intsctysno     = zeros((N+Nactin)*Nsegments,(N+Nactin)*Nsegments);             % an array of 0,1s if the intersection has already happened or no.  
intsctdistside = ones(N*Nsegments,Nsides)*intsctinfty; 

todraw     = 0; % to draw or not to draw
drawtime = 0; % how many times have drawn
cost         = 0; 

cycletime         = zeros(N);
collectcycletime  = zeros(10000); 
cyclenumber       = 0; 
gamma1            = [0,cumsum(gamma(2:end))];   % gamma1 is the array of gammas, where every cell side has the gamma1 angle wrt horizontal, not the previous side. 


timesdrawn = 0; 
sd_array = zeros(10000, 2);

%___ for saving files for plotitng later
state_cell = zeros(N*Nsegments + 1, 6); % Here the 6  elements are : r  cos(theta)    sin(theta) theta   xseed    yseed 

index_saved_figure = 1;
t_save_figure =  [t_end-3 : 0.03 : t_end];



%____________( Main Cycle )_____________________________________________

tic
 while ((time<t_end) ), 
    tindex = tindex+1;
    % _______________( figure out  which site is active )__________
    bigsum  = sum(bigstate);
    [mm,ii] = max( real((cumsum(bigstate./bigsum) - rand)>0)); 
                                                                % the first max gives the active site and what it does.

    ab = mod(ii,2);                                             % if it's alpha or betas: ab<0.5 -> alphas, ab>0.5 -> betas.
    if (ab<0.5), siten = floor(ii)/2; else siten = (ii+1)/2; end % site number 
    

    % --- one growing cycle ---                                 % --- one growing cycle ---    
    if capyesno(siten)>0.5,                                     % if have a cap, then in An
        if ab > 0.5,                                            %    if (alpha  = if have a cap and growing) 
           
            [site, segm, xseed, intsctysno, intsctdist, intsctdistside,capyesno,bigstate,state] ...
         = proga_addonedimer06(site,siten,segm,sidenode,seedside,xseed,gamma,intsctysno,intsctdist,intsctdistside,N,Nsegments,Nactin,L,intsctinfty,capyesno,bigstate,state,thetacrwall,thetacrmt,Pcat);
            
            cost = cost + 1; 
        else                                                    %    else (if betapr = if have a cap and losing a cap)
                                                                %       ['site = ',num2str(siten),' losing a cap, shrinks, rtotal = ', sum(site(siten,:,1))]
           capyesno(siten) = 0; bigstate([2*siten-1, 2*siten]) = state(:,1);    % cap->lost, update bigstate
          [site,segm,capyesno,intsctysno,intsctdist,intsctdistside,bigstate] = ...
          proga_shrinking001(site,siten,segm,Nsegments,capyesno,intsctysno,intsctdist,intsctdistside,bigstate,state,intsctinfty);     
        end                                                     %    end
    else                                                        % else (dont' have a cap--> in Bn)
        if ab > 0.5,                                            %    if   (alphapr  = don't have a cap,   getting a cap and adding a tubulin) 
            capyesno(siten) = 1; bigstate([2*siten-1, 2*siten]) = state(:,2);  % get a cap, update bigstate
            if site(siten,1,1) <0.5,                            %        initialize theta, if starting from zero
               site(siten,1,2:4) = randarchonside(seedside(siten),gamma); % initialize the directions for a new MT site   
               segm(siten) = 1;                                 % do say that we've entered the first segment 
               cycles = cycles + 1; 
            end
            
            [site, segm, xseed, intsctysno, intsctdist, intsctdistside,capyesno,bigstate,state] ...
         = proga_addonedimer06(site,siten,segm,sidenode,seedside,xseed,gamma,intsctysno,intsctdist,intsctdistside,N,Nsegments,Nactin,L,intsctinfty,capyesno,bigstate,state,thetacrwall,thetacrmt,Pcat);
            
           cost = cost + 1; 
        else                                                    %    else ( beta  = don't have a cap and shrinking)
          
           [site,segm,capyesno,intsctysno,intsctdist,intsctdistside,bigstate] = ...
           proga_shrinking001(site,siten,segm,Nsegments,capyesno,intsctysno,intsctdist,intsctdistside,bigstate,state,intsctinfty);
 
        end                                                     %   end 
        
    end                                                         % end 

    % --- end one growing cycle ---
    
    ddt = - log(rand)/bigsum;                                   % time increment
    time = time  + ddt;                                         % adjust the time 
    
    %%%%%%%% Define specific time intervals for collecting data/plotting %%
    %%%================== CALCULATING THE MTSD ============================
    if time > drawtime + dt, 
            drawtime   = drawtime + dt; 
            timesdrawn = timesdrawn + 1;  
            todraw     = 1; 
            avglengthmt(timesdrawn)     = sum(sum(site(:,:,1)))/sum(sum(site(:,:,1)>0));
            % avg angle 
            uu                          =  sum(sum(   site(:,:,1).* asin((site(:,:,3)*V(1,1) -  site(:,:,2)*V(2,1)))     ))./ sum(sum( site(:,:,1)))     ; 
            uustd(timesdrawn)           =  180/pi*sqrt(    sum(sum(   site(:,:,1).* asin((site(:,:,3)*V(1,1) -  site(:,:,2)*V(2,1))).^2  ))./ sum(sum( site(:,:,1)))   ); 
            uustd0 (timesdrawn)         =  180/pi*sqrt(    sum(sum(   site(:,:,1).* (acos(site(:,:,2)) - uu).^2))  ./ sum(sum( site(:,:,1)))               ); 
            avguu(timesdrawn)           = uu*180/pi; 
            drawingtimearray(timesdrawn)=  time;
            mtdensity(timesdrawn)       =  3* sum(sum(site(:,:,1))) / cellarea; 
            fracnonzeroMTs (timesdrawn) = sum(site(:,1,1)>0)/N; 
            
            

            rr            = zeros(N*Nsegments, 1); 
            theta         = zeros(N*Nsegments, 1); 
            nonzero_segmn = 0;  % the number of nonzero segments 

            % remove all the segments with zero length
            for siten = 1:N, 
                for segmn = 1:Nsegments, 
                    if site(siten,segmn,1) > 0.5,  
                        nonzero_segmn  = nonzero_segmn  + 1; 
                        theta(nonzero_segmn) = mod(site(siten,segmn,4),pi) *180/pi; 
                        rr    (nonzero_segmn) =     site(siten,segmn,1); 
                    end
                end
            end
            % bin theta into 1-180 bins AND keep the count of how many there were there... 
            [sorted_theta,ii] = sort(theta(1:nonzero_segmn)); 
            rr(1:nonzero_segmn) = rr(ii); 

            % only need to store the histogram = 180 numbers. every. single. time. 
            theta_hist = zeros(180,1); 
            for jj = 1:nonzero_segmn, 
                kk = ceil(sorted_theta(jj)); 
                if kk ==0, kk = 1; end
                theta_hist(kk) = theta_hist(kk)  + rr(jj); 
            end
            theta_hist_avg = theta_hist_avg * (timesdrawn / (timesdrawn + 1)) + theta_hist * (1/timesdrawn);  % running time-average of the histogram
            
            
            stretch_factor = 0.5;
            bin_size       = 4;
            bin            = (bin_size*(1/stretch_factor) * pi/180); %bin size in radians for Von Mises Fit


            x_direction_radians = [0 : (bin_size*2) : (360 - (bin_size*2))]'*pi/180;  %Takes into account the bin size to build the x values
            number_bins         = max(size(x_direction_radians)); 
            y                   = zeros(number_bins,1);
     
            for ii = 1:number_bins,
                y(ii) = sum(theta_hist_avg((ii-1)*4 + 1 : (ii-1)*4 + 4));
            end 

            kappa = circ_kappa(x_direction_radians, y, bin);
            SD = sqrt(1/kappa) * (180/pi) * stretch_factor;
            mu = circ_mean(x_direction_radians,y);
            mu_degrees = mu * (180/pi) * stretch_factor;
         
              
            sd_array(timesdrawn,:)  = [SD, drawtime];
            
  %%%================== CALCULATING BUNDLE INFO ============================          
  
  % create a matrix which stores the angle for each MT and segm & radii on
  % test line xing for the number of test lines investigated
  xing_testline = zeros(N*Nsegments, 1+lines_measure);
  
  % use biggest test line for normalizing the test lines for frequency
  % plots
  check_max_line = zeros(1,lines_measure);
  for kk = 1:lines_measure,
    check_max_line(kk) = site(N+Nactin+kk,1,1);
  end
  % !!! This is currently not being used for standardizing
  MAX_lines = max(check_max_line);
  angles_all = zeros(1,1);
  radii_all = zeros(1,1);
  
  for jj = 1:lines_measure,
       for siten = 1:N,
         
           if segm(siten) > 0.5;
             for seg = 1:segm(siten),
                   
                   segmn_overall = (siten-1)*Nsegments + seg; 

                   %insert relevant angles into the first column
                   xing_testline(segmn_overall,1) = site(siten,seg,4);

                   x1 = xseed(siten,seg,1);
                   y1 = xseed(siten,seg,2);
                   theta1 = site(siten,seg,4); 

                   x2 = xseed(N+Nactin+jj,1,1);
                   y2 = xseed(N+Nactin+jj,1,2);
                   theta2 = site(N+Nactin+jj,1,4);

                   [r1,r2] = intersectionrs01(x1,y1,theta1,x2,y2,theta2);

                   % r1 will give the critical radii of MT for xing with testline
                   r1 = (r1<0)*intsctinfty + (r1>=0).*r1;  r1 = min(intsctinfty,r1);

                   if r1 <= site(siten,seg,1),
                   % r2 is the distance on the test line
                   % if MT is bigger than critical dist
                   r2 = (r2<0)*intsctinfty + (r2>=0).*r2;  r2 = min(intsctinfty,r2);
                   xing_testline(segmn_overall,1+jj) = r2;
                   % (r2/MAX_lines)*60;
                   else
                   % if MT is smaller than critical dist
                   xing_testline(segmn_overall,1+jj) = 0;
                   end
               end
           else
           end
       
       end
  end
  
  %%%%%%%%%%%%%%%%%%%% BUNDLE CALCULATIONS FOR TEST LINES %%%%%%%%%%%%%%%%%
  % create matrix which stores the bundle info for all test lines plus the
  % overall information!!!
  bundles_info = zeros(3,lines_measure + 1);
  
  
  for mm = 1:lines_measure,
      
      m1 = mm+1;
      
      %select angle info and test line
      nonzero_bund = xing_testline(:,[1 m1]);
      
      % picks out all elements which have two nonzero components, ie.
      % nonzero rows
      nonzero_bund = nonzero_bund(all(nonzero_bund,2),:);
      
      % calculate all active points xing test line
      total_MTS_xing_line = length(nonzero_bund(:,1));
      
      % create matrix for storing radii and frequency of MTS that are in
      % bundle at that point 
      all_bundles = ones(length(nonzero_bund(:,1)),2);
      %%% !!! CHANGE ALL_BUNDLES SO THAT THEY HAVE THE RADII INITIALLY OF
      %%% ALL XING, ONLY SETUP THE COUNTER THINGY...!!!
   
      
      % go through length of nonzero xings with testlines and compare them to all other
      % non zero xings
      for uu = 1:length(nonzero_bund(:,1)),
          
          
          
           % find all indices apart from the current uu
           m = 1:length(nonzero_bund(:,1));
           ss = circshift(m',length(m)-uu);   
           ss = ss(1:end-1)';
     
          for kk = 1:length(ss);     
               % first check if these hit the test line at the same spot
               % and check if angles are multiples of one another
               if abs(nonzero_bund(uu,2)-nonzero_bund(ss(kk),2))==0.3 && pi + pi/180 < abs((nonzero_bund(uu,1)-nonzero_bund(ss(kk),1))) < pi + pi/90 || abs((nonzero_bund(uu,1)-nonzero_bund(ss(kk),1))) < pi/180,
                   
                   % set first column to radii 
                   all_bundles(uu,1) = nonzero_bund(uu,2);
                   % ADD ONE TO NO OF MTS IN BUNDLE
                   all_bundles(uu,2) = all_bundles(uu,2) + 1;
                   
               else       
               end
          end
      end
   
      
    
      kon = find(all_bundles(:,2)>1.5);
      final_bundle_info = all_bundles(kon,2);
      
      % groups bundles into amounts ie. first arranges in ascending order 
      final_bundle_info = sort(final_bundle_info);
      
      % counts how many appear ie. [2 2 2 2 3 3 3] will give unique = [2 3]
      % and histc as [4 3] > take [4 3]./[2 3] to get [2 1] ie 2 bundles of
      % 2 and 1 bundle of size 3
      types_bundles = unique(final_bundle_info);
      number_bundles_vec = histc(final_bundle_info,types_bundles)./types_bundles;
    
      % counts total number of bundles
      number_bundles = sum(number_bundles_vec);
      
      
      % calculate the percentage of bundles
      perc_bundles = (number_bundles/total_MTS_xing_line)*100;
      
      
      % ie. [2 1].*[2 3] = [4 3] sum[4 3]=7/3=2.3 
      average_bundle_size = sum(number_bundles_vec.*types_bundles)/number_bundles;
      
      bundles_info(1,mm) = perc_bundles;
      bundles_info(2,mm) = average_bundle_size;
      bundles_info(3,mm) = number_bundles;
      
  end
  
  
  
% %   %%%%%%%%%%%%%%%%%%%% OVERALL BUNDLE CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%

% First to lower computational cost, check which "overall_segment" indices
% are greater than zero !!! we will be working with a matrix of this
% dimension

bundles = 0;
for i = 1:N,
    if segm(i) > 0.5,
         jj = [1:segm(i)];
         segmn_overall = (i-1)*Nsegments + jj;
            
            if bundles == 0,
                   bundles = segmn_overall;
            else
                   bundles = [bundles segmn_overall];
            end
    end
end

BUNDLES_ALL = zeros(length(bundles),length(bundles));
BUNDLES_ALL(:,1) = bundles';

% index_brightness = proga_index_brightness(N,xseed,segm,site,half_width);
% index_brightness_segment = sum(index_brightness,2); 

for i = 1:N,
    if segm(i) > 0.5,
        for j = 1:segm(i),
            segmn_overall = (i-1)*Nsegments + j;
            a = find(BUNDLES_ALL(:,1) == segmn_overall);
            
            % find all other MTs that have to be checked against
            b = 1:N;
            cc = circshift(b',length(b)-i);
            cc = cc(1:end-1);
            
            for i2 = 1:N-1,
                if segm(cc(i2)) > 0.5,
                    
                    counter = 0;
                    
                    for j2 = 1:segm(cc(i2)),
                        segmn_overall2 = (i2-1)*Nsegments + j2;
                        b = find(BUNDLES_ALL(:,1) == segmn_overall2);
                        
%                         % check if the two segments have zipped using
%                         % brightness code
%                         if index_brightness_segment(segmn_overall) > 0.5 && index_brightness_segment(segmn_overall2) > 0.5,
                            
                            % if they have zipped check if their angles are
                            % the same ie. if they are in the same bundle
                            if abs(site(i,j,4)-site(i2,j2,4)) == pi || abs(site(i,j,4)-site(i2,j2,4)) == 0,
                                
                                
                                                 if j-1 > 0.5,

                                                 x1 = xseed(i,j-1,1); y1 = xseed(i,j-1,2);  theta1  = site(i,j-1,4);  
                                                 x2 = xseed(i2,j2,1); y2 = xseed(i2,j2,2);  theta2  = site(i2,j2,4);  

                                                 [r1,r2] = intersectionrs01(x1,y1,theta1,x2,y2,theta2);
                                                 r1 = (r1<0)*intsctinfty + (r1>=0).*r1;  r1 = min(intsctinfty,r1);
                                                 r2 = (r2<0)*intsctinfty + (r2>=0).*r2;  r2 = min(intsctinfty,r2);

                                                     if site(i,j-1,1) + 1 > r1,
                                                         counter = counter + 1;
                                                     else
                                                         counter = counter;
                                                     end


                                                 elseif j2-1 > 0.5,

                                                 x1 = xseed(i,j,1); y1 = xseed(i,j,2);  theta1  = site(i,j,4);  
                                                 x2 = xseed(i2,j2-1,1); y2 = xseed(i2,j2-1,2);  theta2  = site(i2,j2-1,4);  

                                                 [r1,r2] = intersectionrs01(x1,y1,theta1,x2,y2,theta2);
                                                 r1 = (r1<0)*intsctinfty + (r1>=0).*r1;  r1 = min(intsctinfty,r1);
                                                 r2 = (r2<0)*intsctinfty + (r2>=0).*r2;  r2 = min(intsctinfty,r2);

                                                     if site(i2,j2-1,1) + 1 > r2,
                                                            counter = counter + 1;
                                                     else
                                                         counter = counter;
                                                     end

                                                         if counter > 0.5,
                                                             BUNDLES_ALL(a,b) = segmn_overall2;
                                                         end

                                                 else
                                                      BUNDLES_ALL(a,b) = 0;
                                                 end

                            else
                                BUNDLES_ALL(a,b) = 0;
                            end
                    end
                end
            end
        end
    end
end

BUNDLES_ALL = sort(BUNDLES_ALL,2);
bundles2 = 0;

for kk = 1:length(BUNDLES_ALL(:,1)),
    
    if length(nonzeros(BUNDLES_ALL(kk,:))) > 1.5,
        if bundles2 ==0,
            
            bundles2 = BUNDLES_ALL(kk,:);
        else
            bundles2 = [bundles2 ; BUNDLES_ALL(kk,:)];
        end
    end
    
end

number_of_bundles_all = length(unique(bundles2(:,end)));

no_all_nonzero_bundles = length(bundles2(:,1));
percentage_all_bundles = (no_all_nonzero_bundles/nonzero_segmn)*100;


% 
%                          counter1 = 0;
%                          counter2 = 0;
%                          
%                          for seg2 = 1:segm(cc(siten2)),
% 
%                              
%                                   if abs(site(siten,seg,4)-site(cc(siten2),seg2,4)) == pi || abs(site(siten,seg,4)-site(cc(siten2),seg2,4)) == 0,
% 
%                                       if seg2 > 1.5,
%              x1 = xseed(siten,seg,1); y1 = xseed(siten,seg,2);  theta1  = site(siten,seg,4);  
%              x2 = xseed(cc(siten2),seg2-1,1); y2 = xseed(cc(siten2),seg2-1,2);  theta2  = site(cc(siten2),seg2-1,4);  
% 
%              [r1,r2] = intersectionrs01(x1,y1,theta1,x2,y2,theta2);
%              r1 = (r1<0)*intsctinfty + (r1>=0).*r1;  r1 = min(intsctinfty,r1);
%              r2 = (r2<0)*intsctinfty + (r2>=0).*r2;  r2 = min(intsctinfty,r2);
%                           if r2 <= 1 + site(cc(siten2),seg2-1,1)  && r1 <= site(siten,seg,1),
%                               counter1 = counter1 + 1;
%                           end
%                                       end
% 
%                                       if seg > 1.5,
%              x1 = xseed(siten,seg-1,1); y1 = xseed(siten,seg-1,2);  theta1  = site(siten,seg-1,4);  
%              x2 = xseed(cc(siten2),seg2,1); y2 = xseed(cc(siten2),seg2,2);  theta2  = site(cc(siten2),seg2,4);  
% 
%              [r1,r2] = intersectionrs01(x1,y1,theta1,x2,y2,theta2);
%              r1 = (r1<0)*intsctinfty + (r1>=0).*r1;  r1 = min(intsctinfty,r1);
%              r2 = (r2<0)*intsctinfty + (r2>=0).*r2;  r2 = min(intsctinfty,r2);
%                           if r2 <= site(cc(siten2),seg2,1)  && r1 <= 1 + site(siten,seg-1,1) ,
%                               counter2 = counter2 + 1;
%                           end

       bundles_info(1,end) = percentage_all_bundles;
%       bundles_info(2,end) = average_bundle_size_total;
       bundles_info(3,end) = number_of_bundles_all;
      
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %Bin into 1-60 for the frequency plot, middle line only !!! =============
  bins_line = 0:60;
  
  frequency_line = nonzeros(xing_testline(:,1+(lines_measure+1)/2));
  %since seeds are on the rhs flip this around !!!
  frequency_line =  site(N+Nactin+((lines_measure+1)/2),1,1) - frequency_line;
  freq_line = histc(frequency_line,bins_line);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
            %______________ saving the cell state once in a while________
   
             if and(time > t_save_figure(index_saved_figure), time < t_save_figure(index_saved_figure) + dt),
                index_saved_figure = index_saved_figure  + 1; 
                
                for siten = 1:N, 
                    for segmn = 1:Nsegments, 
                        state_cell(siten*(Nsegments - 1) + segmn, : ) = ...
                           [site(siten,segmn,1) ...
                           site(siten,segmn,2) ...
                           site(siten,segmn,3) ...
                           site(siten,segmn,4) ...
                           xseed(siten,segmn,1)...
                           xseed(siten,segmn,2)];
                    end
                end
                state_cell(N*Nsegments + 1,:) = [N Nsegments time L(1) 0 0];
                
                filename_state_cell = [dir_name,'state_cell_',num2str(ecc*100),'_t_',num2str(index_saved_figure),'_testnumber_',num2str(testnumber),'_Actin_',num2str(Nactin),'_ap4bp1a1000b3.5aN220L120_2.csv'];
                csvwrite(filename_state_cell,          state_cell);
                disp(['____ saved figure to',filename_state_cell]);

                filename_bundles = [dir_name,'bundles_',num2str(ecc*100),'_t_',num2str(index_saved_figure),'_testnumber_',num2str(testnumber),'_Actin_',num2str(Nactin),'_ap4bp1a1000b3.5aN220L120_2.csv'];
                csvwrite(filename_bundles,  bundles_info);
                
             end

     end
                   
          
     if draw0screen1 == 0, 
         if and( todraw > 0.5, time>timetodraw),

             
   %%%%% plotting with actin uses the site information to plot the Actin seeds          
               proga_drawing08_actin(theta_hist_avg, xseed,segm,sidenode, L, N,Nactin,site,time,avglengthmt,uu,SD,uustd0,avguu,...
                   drawingtimearray,V,ecc,avglengthmttheory,xseed1,mtdensity,fracnonzeroMTs,timesdrawn,half_width,dir_name,lines_measure,MAX_lines,bins_line,freq_line,bundles_info);
               todraw = 0; 
         end 
     else
         if and( todraw > 0.5, time>timetodraw),
            disp(['t = ', num2str(time),', SD = ', num2str( SD ),', frac nonzero mts  = ', num2str(fracnonzeroMTs (timesdrawn)), ...
              ', avg(length MT) = ', num2str(   avglengthmt(timesdrawn))]);
            todraw = 0; 
        end
     end
     
  
     
 end

   
toc        

site_for_juan = zeros(N*Nsegments,2); 

for siten = 1:N, 
    for segmn = 1:Nsegments, 
        osegmn        = (siten-1)*Nsegments + segmn;
        site_for_juan(osegmn,1) = site(siten,segmn,4);
        site_for_juan(osegmn,2) = site(siten,segmn,1);
    end
end

qwerqwer = theta_hist_avg; 



filenamesused_theta_hist_avg = [dir_name,'angles_lengths_ecc_',num2str(ecc*100),'testnumber_',num2str(testnumber),'Actin_',num2str(Nactin),'_ap4bp1a1000b3.5aN220L60_1.csv']
filenamesused_mtsd = [dir_name,'mtsd_time_evolution',num2str(ecc*100),'testnumber_',num2str(testnumber),'Actin_',num2str(Nactin),'_ap4bp1a1000b3.5aN220L60_1.csv']


 csvwrite(filenamesused_theta_hist_avg,theta_hist_avg);
 csvwrite(filenamesused_mtsd,          sd_array);
 
   %______________ saving the cell state once in a while________
   
                index_saved_figure = index_saved_figure  + 1; 
                
                for siten = 1:N, 
                    for segmn = 1:Nsegments, 
                        state_cell(siten*(Nsegments - 1) + segmn, : ) = ...
                           [site(siten,segmn,1) ...
                           site(siten,segmn,2) ...
                           site(siten,segmn,3) ...
                           site(siten,segmn,4) ...
                           xseed(siten,segmn,1)...
                           xseed(siten,segmn,2)];
                    end
                end
                state_cell(N*Nsegments + 1,:) = [N Nsegments time L(1) 0 0];
                
                filename_state_cell = [dir_name,'state_cell_',num2str(ecc*100),'_t_',num2str(index_saved_figure),'testnumber_',num2str(testnumber),'Actin_',num2str(Nactin),'_ap4bp1a1000b3.5aN220L60_1.csv'];
                csvwrite(filename_state_cell,          state_cell);

   
                filename_angles_all = [dir_name,'ANGLES_ALL_',num2str(ecc*100),'_testnumber_',num2str(testnumber),'_Actin_',num2str(Nactin),'_ap4bp1a1000b3.5aN220L120_2.csv'];
                csvwrite(filename_angles_all, angles_all_overall);
 
end
  