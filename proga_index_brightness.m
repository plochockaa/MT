% the function proga_index_brightness.m gives the brightness of each
% segment of the MT. 
% 
% The brightness index is the number of >>> MTs <<< that are "nearby", 
% but are not crossing any given segment. 
%
% The code for each segment, 
% -- checks if the angle with any other MT segment <30, 
% -- if it's <30, then we check if the segments are nearby or not (wrt "fat" MT of width 3 dimers
% -- 
%
% INPUT 
%     N, Nsegments = overall # of MTs and MT segments
%     site  = the array with all the info about the segments, angles,etc
%     xseed = the array with the info about the MT seeds
%     segm  =  
%
% OUTPUT
%
%
%


function   index_bightness = proga_index_brightness(xseed,segm,site,half_width);


    qwer = size(site); 
    Nsegments = qwer(2); 
    N         = qwer(1);
    
    gamma_critical = 1/180*pi;

    index_bightness = zeros(N*Nsegments,N); 


    for siten1 = 1:N-1, 
    if segm(siten1) >0.5,     
    for segmn1 = 1:segm(siten1), 
            
            osegmn1 = (siten1-1)*Nsegments + segmn1;       
            
            for siten2 = siten1+1:N,
            if segm(siten1) >0.5,     
            for segmn2 = 1:segm(siten2), 
            if or ( or ( abs(site(siten1,segmn1,4)      - site(siten2,segmn2,4)) < gamma_critical,  ...
                         abs(site(siten1,segmn1,4) + pi - site(siten2,segmn2,4)) < gamma_critical ), ...
                         abs(site(siten1,segmn1,4) - pi - site(siten2,segmn2,4)) < gamma_critical ), 
                
               
                %index_inside = nearby_or_not(seed_coords_1, r1, cos_theta1, sin_theta1,  seed_coords_2, r2, cos_theta2, sin_theta2, half_width);
                index_inside = nearby_or_not(...
                   [xseed(siten1,segmn1,1); xseed(siten1, segmn1,2)], ...   % seed_coords_1, 
                    site (siten1,segmn1,1), ...                             % r1, 
                    site (siten1,segmn1,2), ...                             % cos_theta1, 
                    site (siten1,segmn1,3), ...                             % sin_theta1,  
                   [xseed(siten2,segmn2,1); xseed(siten2, segmn2,2)], ...   % seed_coords_2, 
                    site (siten2,segmn2,1), ...                             % r2, 
                    site (siten2,segmn2,2), ...                             % cos_theta2, 
                    site (siten2,segmn2,3), ...                             % sin_theta2, 
                    half_width);
                
                 osegmn2 = (siten2-1)*Nsegments + segmn2;  
                 
                 index_bightness(osegmn1,siten2) =  index_bightness(osegmn1,siten2) + index_inside ;
                 index_bightness(osegmn2,siten1) =  index_bightness(osegmn2,siten1) + index_inside ; 
                     
        
            end 
            end
            end
            end
    end
    end
    end


end

