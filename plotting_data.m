
 clear all; close all; clc;
 
 
 %ecc_array = [0.7 0.75 0.8 0.85 0.9 0.95];
 
 mtsd_array = ecc_array*0;
 figure(2); clf; 
 
    stretch_factor = 0.5;
    bin_size       = 4;
    bin            = (bin_size*(1/stretch_factor) * pi/180); %bin size in radians for Von Mises Fit


    x_direction_radians = [0 : (bin_size*2) : (360 - (bin_size*2))]'*pi/180;  %Takes into account the bin size to build the x values
    number_bins         = max(size(x_direction_radians)); 
    y                   = zeros(number_bins,1);

 

 for jj = 1:length(ecc_array), 
    ecc = ecc_array(jj); 
    
    % read the avg of angle histogram and the array of MTSD versus time (time is the 2nd column)
    filenamesused = ['angles_lengths_ecc_',num2str(ecc*100),'_ap4bp1a1000b3.5aN220.csv'];
  % filenamesused_mtsd = ['mtsd_time_evolution',num2str(ecc*100),'_ap4bp1a1000b3.5aN220.csv'];
    
    
    
% for the lentgh 60, these are the file names: 
filenamesused_theta_hist_avg = ['angles_lengths_ecc_',num2str(ecc*100),'_ap4bp1a1000b3.5aN220L60_1.csv'];
filenamesused_mtsd           = ['mtsd_time_evolution',num2str(ecc*100),'_ap4bp1a1000b3.5aN220L60_1.csv'];



    theta_hist_avg = csvread(filenamesused_theta_hist_avg);
    sd_array       = csvread(filenamesused_mtsd);    
 
    figure(1); 
        subplot(3,3,jj); 
            a = sd_array;  
            a(find(a(:,1)==0), :) = []; % delete zeros
            plot(a(:,2), a(:,1),'LineWidth',1.5); hold on; 
            mean_mtsd = mean(a(floor(length(a)*0.75):end,1)); 
            plot([0:max(a(:,2))],mean_mtsd + [0:max(a(:,2))]*0,'r','LineWidth',1.5);
            axis([0 10 25 70])
            title(['ecc = ', num2str(ecc),',  MTSD -> ', num2str(mean_mtsd)],'FontSize', 16)
            ylabel('MTSD','FontSize', 16); xlabel('time','FontSize', 16); 
            grid on; 
            
            
            mtsd_array(jj)= mean_mtsd; 
     
            
 end
    figure(4); plot(ecc_array, mtsd_array,'o-','Markersize',5,'LineWidth',2); grid on; 
    title('MTSD(ecc)','FontSize', 15) 
    xlabel('ecc','FontSize', 15) 
    ylabel('MTSD','FontSize', 15)
 
 for jj = 1:length(ecc_array)-1, 
    ecc = ecc_array(jj); 
    
    % read the avg of angle histogram and the array of MTSD versus time (time is the 2nd column)
    filenamesused = ['angles_lengths_ecc_',num2str(ecc*100),'_ap4bp1a1000b3.5aN220.csv'];
    theta_hist_avg = csvread(filenamesused);
            
    figure(2); 
        subplot(2,3,jj); 
            for ii = 1:number_bins,
                y(ii) = sum(theta_hist_avg((ii-1)*4 + 1 : (ii-1)*4 + 4));
            end 

            kappa      = circ_kappa(x_direction_radians, y, bin);
            SD         = sqrt(1/kappa) * (180/pi) * stretch_factor;
            mu         = circ_mean(x_direction_radians,y);
            mu_degrees = mu * (180/pi) * stretch_factor;

            xx = [1:180]; 

            fitted_distribution = exp(kappa * cos(xx*pi/180/stretch_factor - mu));
            fitted_distribution = fitted_distribution / norm(fitted_distribution); 

            hold on; 
            plot(xx,theta_hist_avg/norm(theta_hist_avg),'b'); 
            plot(xx,fitted_distribution,'m','LineWidth',2); 

            title(['K = ', num2str(kappa),',  SD = ',num2str(SD),',  mu = ', num2str(mu_degrees)],'FontSize',15)
            axis([0 180 0 1/4])

        
 end
 
 
 

 
 
 %%%%%%%%%%%%%%%% last part, plotting all the cells %%%%%%%%%%%%%%%%
 
 
    
    % main plotting cycle 
    for index_saved_figure = 82, %95 ,
        
    for jj = 1:length(ecc_array), 
    ecc = ecc_array(jj); 
    
    % build a cell 
    a   = 1; 
    b   = a./sqrt(1 - ecc^2); 
    L     = [a b  a b]; 
    gamma = [0   80 110 0  ]/180*pi;
  
  
    Nsides   = length(L); 
   
         
    
  
    filename_state_cell = ['state_cell_',num2str(ecc*100),'_t_',num2str(index_saved_figure),'_ap4bp1a1000b3.5aN220L120_2.csv'];
    state_cell = csvread(filename_state_cell);

    N          =  state_cell(end,1);
    Nsegments  =  state_cell(end,2);
    time       =  state_cell(end,3);
    L0          =  state_cell(end,4);
    
      L         = L*L0; 
    
    
     xseed1    = zeros(N,2);   % N = number of MTs, xseed = seed positions of the N MTs
    seedsdie = []; 
    sidenode = []; 
   
         

    [xseed1,seedside,L,gamma,sidenode] = nodes_initiation(N,Nsegments,L,gamma);
    
     xseed    = zeros(N,Nsegments,2);
    
    
    
   
    figure_number = index_saved_figure; 
    figure(figure_number); 
    subplot(2,4,jj)
       hold on; 
        plot([sidenode(:,1);  sidenode(1,1)], [sidenode(:,2);  sidenode(1,2)],'Color',[1 0 1],'LineWidth',3, 'LineSmoothing','on'); % plot the cell sides 
       % whitebg(1,'k')
       
        
        % find the siten,segmn array
      
         
         xseed  = zeros(N,Nsegments,2);
         site   = zeros(N,Nsegments,4);
         
%          full_index = [1:N*Nsegments]; 
          number_of_neighbors = zeros(N*Nsegments, 1)+1; 
         
         for osegmn = 1:N*Nsegments, 
            segmn = mod(osegmn, Nsegments);
            siten = (osegmn - segmn)/Nsegments + 1;
             if segmn == 0, segmn = Nsegments; siten  = siten  - 1; end      
            
            site(siten,segmn,  1) = state_cell(osegmn,1);
            site(siten,segmn,  2) = state_cell(osegmn,2);
            site(siten,segmn,  3) = state_cell(osegmn,3);
            site(siten,segmn,  4) = state_cell(osegmn,4);
       
            xseed(siten,segmn,1) = state_cell(osegmn,5); 
            xseed(siten,segmn,2) = state_cell(osegmn,6);
            
%             for kk = min(osegmn+1,Nsegments-1):Nsegments-1,
%                 
%                 if  state_cell(kk,1) > 0.5, 
%                 if ((state_cell(osegmn,5) -  state_cell(kk,5))^2 + (state_cell(osegmn,6) - state_cell(kk,6))^2 < 9), 
%                     number_of_neighbors(osegmn) = number_of_neighbors(osegmn) + 1;
%                     number_of_neighbors(kk) = number_of_neighbors(kk) + 1;
%                 end
%                 end 
%             end
            
         end
%          if sum(number_of_neighbors) ==0, 
%              'FAAAK'
%             pause
%          end
%          
%          number_of_neighbors = number_of_neighbors/max(1, max(number_of_neighbors));
         
        beh = max(L)/10; 
        xlim([min(xseed(:,1,1))  - beh, max(xseed(:,1,1)) +  beh]);
        ylim([min(xseed(:,1,2))  - beh, max(xseed(:,1,2)) +  beh]);

        
        alpha_plotting  = 0.4; 
        
        for siten =1:N, 
            for segmn = 1:Nsegments, 
                if site(siten,segmn,1) >0.5, 
                       
                   
                    
                    
                        x = site(siten,segmn,2)*site(siten,segmn,1);    %correct, since the r = uvelichenie :) and x and y are the initial values on a unit circle. omg. i was so smart. 
                        y = site(siten,segmn,3)*site(siten,segmn,1);
                        
                       colorn = alpha_plotting  + (1-alpha_plotting) *  number_of_neighbors((siten - 1)*Nsegments + segmn); 
                        
                        plot([xseed(siten,segmn,1) , xseed(siten,segmn,1) + x]',...
                             [xseed(siten,segmn,2) , xseed(siten,segmn,2) + y]','Color',[0 colorn .1],'LineWidth',1, 'LineSmoothing','off'); hold on; 
               
                end 
            end
        end
        
        
        %__________( begin: find mtsd )____________________
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
            theta_hist_avg = theta_hist; % no averaged histogram this time cause need just one instance in time.  
            
            
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
        %__________( end  : find mtsd )____________________
       
        axis equal
       % whitebg(1,'k')
      set(gca,'Color',[0 0 0]);
   
        title(['ecc=',num2str(ecc), ' mtsd=',num2str(SD) ],'FontSize',15)
    
    
    
    
   
    end 

        
 end
 
 
 
 

 