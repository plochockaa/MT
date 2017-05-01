function                proga_drawing08_actin(theta_hist_avg, xseed,segm,sidenode, L, N,Nactin,site,time,avglengthmt,SD,uustd0,avguu,...
                   drawingtimearray,V,ecc,avglengthmttheory,xseed1,mtdensity,fracnonzeroMTs,timesdrawn,half_width,dir_name);

    qwer =size(site(1:N,:,:)); 
    Nsegments = qwer(2); 
               

    index_bightness = proga_index_brightness(N,xseed,segm,site,half_width);
    index_bightness_segment = sum(index_bightness,2); 

    color_factor = 0.8/6; % max(1,max(index_bightness_segment)); 

    figure(1); 
   clf
   set(gcf, 'Position', [200   200   500   500]);
   %  >>> NOTE <<< 
   % use  get(gcf,'position') to find where the current screen is -- set it manually
  % movegui(figure(1),'northwest');
   

   % subplot(3,4,[1,2,5,6,9,10]);
      %  plot([xseed(:,1,1);xseed(1,1,1)],[xseed(:,1,2);xseed(1,1,2)],'ro'); % plot the cell itself with red circles being the seeds. 
      

          
 
for a=N+1:Nactin+N;
    xactin=site(a,1,1)*cos(site(a,1,4));
    yactin=site(a,1,1)*sin(site(a,1,4)); 
          
    plot([0 xactin]+xseed(a,1,1),[0 yactin]+xseed(a,1,2),'Color',[.4 0.5 1],'LineWidth',.5);
   
hold on
end


  plot([sidenode(:,1);  sidenode(1,1)], [sidenode(:,2);  sidenode(1,2)],'Color',[1 0 1],'LineWidth',3); % plot the cell sides 
        hold on; 


        
        beh = max(L)/10; 
        xlim([min(xseed(:,1,1))  - beh, max(xseed(:,1,1)) +  beh]);
        ylim([min(xseed(:,1,2))  - beh, max(xseed(:,1,2)) +  beh]);

        for i =1:N, 
            if segm(i) > 0.5, 
                for j = 1:segm(i),
                    x = site(i,j,2)*site(i,j,1);    %correct, since the r = uvelichenie :) and x and y are the initial values on a unit circle. omg. i was so smart. 
                    y = site(i,j,3)*site(i,j,1);
%                    plot([xseed(i,j,1) , xseed(i,j,1) + x]',...
%                         [xseed(i,j,2) , xseed(i,j,2) + y]','Color','b','LineWidth',1.2, 'LineSmoothing','on');
%                      plot([xseed(i,j,1) , xseed(i,j,1) + x]',...
%                          [xseed(i,j,2) , xseed(i,j,2) + y]','Color',[0 1 .1],'LineWidth',1, 'LineSmoothing','on');
                     % argh =   color_factor*min(6,index_bightness_segment((i-1)*Nsegments + j)) ; 
                     % argh_color = [0  1 - argh .1 ]; 

                     argh = color_factor *2; 
                     argh_color = [0  1 - argh .1 ]; 
                     
                     % argh_linewidth = 1.5 + min(7,index_bightness_segment((i-1)*Nsegments + j))/7*4;

                     argh_linewidth = 1.5; 
                     
                      plot([xseed(i,j,1) , xseed(i,j,1) + x]',...
                          [xseed(i,j,2) , xseed(i,j,2) + y]','Color',argh_color,'LineWidth',argh_linewidth);
                end
               % plot(xseed(i,segm(i),1) + site(i,segm(i),1)*site(i,segm(i),2), ... 
               %      xseed(i,segm(i),2)  + site(i,segm(i),1)*site(i,segm(i),3),'.r','Markersize',5); % coordinate of the last MT segment + r of that segment 
                 
            else
              %  plot(xseed(i,1,1),xseed(i,1,2),'*g'); % coordinate of the last MT segment + r of that segment 
            end
        end
       
        cc = mean(xseed1); 
        LL = max(max(L)); 
       %  plot(cc(1)  + [-L/2:1:L/2]*V(1,1), cc(2) + [-L/2:1:L/2]*V(2,1),'-m','LineWidth',5);
     %  for ii = -200:1:200,  plot(cc(1)  + ii*V(1,1), cc(2) + ii*V(2,1),'-m','LineWidth',5);end

       
      %  title(['ECC = ',num2str(ecc),', avg angle = ', num2str(180/pi*uu),',  time = ',num2str(time)],'FontSize',15);             
        title(['ECC = ',num2str(ecc),',  time = ',num2str(time)],'FontSize',15);             
        hold on;    axis equal; 
        
        
% movegui(gcf,[0 0]);
 %set(gcf, 'Position', [0 0 500 500 ]);
 text1 =  ['MTSD = ', num2str(SD)];
 % -- used for 0.9 -- text(80,10, text1,'Color','k','FontSize',20,'FontWeight','bold','HorizontalAlignment','left');
 % -- used for 0.85 -- text(70,10, text1,'Color','k','FontSize',20,'FontWeight','bold','HorizontalAlignment','left');
 text(63,0.98, text1,'Color','k','FontSize',15,'FontWeight','bold','HorizontalAlignment','left');
 set(gca,'color','white')
 
%  argh_color = zeros(7,3); 
%  color_factor = 0.8/6;
%   for jj = 1:7, 
%      argh = 1 -  color_factor*(jj-1);
%      argh_color(jj,:) = [0 argh .1]; 
%  end
%  colormap(argh_color);
%  colorbar
%  caxis([0 7]);
%  colorbar('OuterPosition',[0.21 0.65 0.03 0.25],'FontSize',12,'FontWeight','bold','ytick',0.5:1:6.5,'yticklabel',0:1:6,'TickDir','out','TickLength',[0 0 ])
%  
 
% pos=get(gca,'pos');
% set(gca,'pos',[pos(1)*0.95  pos(2) pos(3)*0.95  pos(4)]);
 
 %pos=get(gca,'pos');
 % hc=colorbar('location','east','position',[pos(1) pos(2)+pos(4) pos(3) 0.03]);
 %set(hc,'xaxisloc','east');
  
 % set(gcf,'Color','black'); %whitebg(1,'k'); 
 %>>>>> SAVING THE FIGURE <<<< 
 
 pause(.1)

  figure_filename = [dir_name,'movie_mt_white_time_thick_lines_smoothed_brightness',num2str(10000+timesdrawn),'.png'];
  f = getframe;  imwrite(f.cdata, figure_filename);
 
 




%   
%         subplot(3,4,[3,4]);
%             plot(drawingtimearray(1:timesdrawn), avglengthmt(1:timesdrawn),'-.'); hold on; 
%             plot(drawingtimearray(1:timesdrawn), mean(avglengthmt(1:timesdrawn)) + avglengthmt(1:timesdrawn).*0,'-k');
%             title(['avg MT length, in theory = ',num2str(avglengthmttheory),', blue - experiment'],'FontSize',15); 
%             whitebg(1,'k')
% 
%         subplot(3,4,[7,8]);
%             stretch_factor = 0.5;
%             bin_size       = 4;
%             bin            = (bin_size*(1/stretch_factor) * pi/180); %bin size in radians for Von Mises Fit
% 
% 
%             x_direction_radians = [0 : (bin_size*2) : (360 - (bin_size*2))]'*pi/180;  %Takes into account the bin size to build the x values
%             number_bins         = max(size(x_direction_radians)); 
%             y                   = zeros(number_bins,1);
%      
%             for ii = 1:number_bins,
%                 y(ii) = sum(theta_hist_avg((ii-1)*4 + 1 : (ii-1)*4 + 4));
%             end 
% 
%             kappa = circ_kappa(x_direction_radians, y, bin);
%             SD = sqrt(1/kappa) * (180/pi) * stretch_factor;
%             mu = circ_mean(x_direction_radians,y);
%             mu_degrees = mu * (180/pi) * stretch_factor;
%          
%                 hold on; 
%                 xx = [1:180]; 
% 
%                 plot(xx,theta_hist_avg/sum(theta_hist_avg),'r'); 
%                 plot(xx,exp(kappa * cos(xx*pi/180/stretch_factor - mu))/(2*pi*besseli(0,kappa))/sum(exp(kappa * cos(xx*pi/180 - mu))/(2*pi*besseli(0,kappa))),'m','LineWidth',2);
%                 title(['kappa = ', num2str(kappa),',  SD = ',num2str(SD),',  mu = ', num2str(mu_degrees)],'FontSize',15)
% 
% 
% 
% 
% 
%         subplot(3,4,[11,12]);
%          plot(drawingtimearray(1:timesdrawn), mtdensity(1:timesdrawn),'-m'); hold on; 
%          plot(drawingtimearray(1:timesdrawn), mtdensity(1:timesdrawn).*0,'-r');  
%          plot(drawingtimearray(1:timesdrawn),fracnonzeroMTs(1:timesdrawn),'-.r');
%          legend('mtdensity','zero','fraction nonzero MTs','Location','NorthWest');
%          title('MT density'); 
%          whitebg(1,'k')
%          
        
        pause(0.01); 
        
        hold off; 

        %             rect = get(gcf,'Position');
        %             rect(1:2) = [0 0];
        %MM(counter) = getframe(gcf,rect);

      %  counter =  counter + 1; 
         
end 