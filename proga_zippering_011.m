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
function qwerqwer  = proga_zippering_011...
               (timetodraw, dt, L, gamma, Nsides, xseed, seedside, sidenode, ...
                N,Nsegments,alpha, alphapr, beta, betapr, t_end, drawfinal,drawall, casestudy,...
                thetacrwall,thetacrmt,ecc,avglengthmttheory,xseed1,Pcat,draw0screen1,half_width,dir_name);
            


cycles   = 0;                            % number of cycles that has started
state    = [alphapr, alpha,  alphapr;    % the first column is all the possible states when the cap is off
               beta, betapr, 0     ];    % the second column is all the possible states when the cap is on. 
                                         % the third column is the rate of exiting the zero length situation                     
site     = zeros(N,Nsegments,4);         % initialize all the sites have length zero
segm     = zeros(N,1); 
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
intsctdist     = ones (N*Nsegments,N*Nsegments)*intsctinfty;  % matrix of intersections A(i,j) = if (i) intersects (j), distance for (i) to the intersection. 
intsctysno     = zeros(N*Nsegments,N*Nsegments);             % an array of 0,1s if the intersection has already happened or no.  
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
         = proga_addonedimer06(site,siten,segm,sidenode,seedside,xseed,gamma,intsctysno,intsctdist,intsctdistside,N,Nsegments,L,intsctinfty,capyesno,bigstate,state,thetacrwall,thetacrmt,Pcat);
            
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
         = proga_addonedimer06(site,siten,segm,sidenode,seedside,xseed,gamma,intsctysno,intsctdist,intsctdistside,N,Nsegments,L,intsctinfty,capyesno,bigstate,state,thetacrwall,thetacrmt,Pcat);
            
           cost = cost + 1; 
        else                                                    %    else ( beta  = don't have a cap and shrinking)
          
           [site,segm,capyesno,intsctysno,intsctdist,intsctdistside,bigstate] = ...
           proga_shrinking001(site,siten,segm,Nsegments,capyesno,intsctysno,intsctdist,intsctdistside,bigstate,state,intsctinfty);
 
        end                                                     %   end 
        
    end                                                         % end 

    % --- end one growing cycle ---
    
    ddt = - log(rand)/bigsum;                                   % time increment
    time = time  + ddt;                                         % adjust the time 
    
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
                
                filename_state_cell = [dir_name,'state_cell_',num2str(ecc*100),'_t_',num2str(index_saved_figure),'_ap4bp1a1000b3.5aN220L120_2.csv'];
                csvwrite(filename_state_cell,          state_cell);
                disp(['____ saved figure to',filename_state_cell]);

             end

     end
                   
          
     if draw0screen1 == 0, 
         if and( todraw > 0.5, time>timetodraw),

               proga_drawing08_brightness(theta_hist_avg, xseed,segm,sidenode, L, N,site,time,avglengthmt,uu,SD,uustd0,avguu,...
                   drawingtimearray,V,ecc,avglengthmttheory,xseed1,mtdensity,fracnonzeroMTs,timesdrawn,half_width);
     %           disp(['t = ', num2str(time),', SD = ', num2str( SD),', frac nonzero mts  = ', num2str(fracnonzeroMTs (timesdrawn)), ...
     %         ', avg(length MT) = ', num2str(   avglengthmt(timesdrawn))]);
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



filenamesused_theta_hist_avg = [dir_name,'angles_lengths_ecc_',num2str(ecc*100),'_ap4bp1a1000b3.5aN220L60_1.csv']
filenamesused_mtsd = [dir_name,'mtsd_time_evolution',num2str(ecc*100),'_ap4bp1a1000b3.5aN220L60_1.csv']


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
                
                filename_state_cell = [dir_name,'state_cell_',num2str(ecc*100),'_t_',num2str(index_saved_figure),'_ap4bp1a1000b3.5aN220L60_1.csv'];
                csvwrite(filename_state_cell,          state_cell);

   

 
end
  