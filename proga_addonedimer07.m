% function [site, bigstate] = proga_addonedimer02(site,siten,segmn,0,site(siten,segmn,4)); 
%
% this function grows one unit from an existing point. 
% rold is how long the segment was before adding the new unit.
% if rold = 0, that means that the array is being initiated. 
% if rold >0,  the MT is already established. 
%
% DIFFERENCE FROM THE OTHER CODES:  
%   >>> the code version 05 has MTs COLLAPSE upon touching the wall <<< 
% 
% INPUT VARIABLES: 
%  site -- all the sites
%  siten -- the number of MT
%  segmn -- the number of the segment 
%  segm(N,1) = for every MT gives the number of nonzero segments the MT has
%  rold = the number of dimers in the segment 
%  dirdata = [x,y,theta] of the direction of the old segment in which the dimer is *expected* to grow  
%            (unless there is zipping along walls or mts).
%  intsctysno = (N*Nsegmn, N*Nsegmn)
%  intsctdist = (N*Nsegmn, N*Nsegmn)
%  intsctdistside = (N*Nsegmn,Nsides)
%
% OUTPUT VARIABLES: 
%  site -- update version of site, where the MT is either grown or not. 
%  tolosecap = 1 (if should lose a cap and change bigstate)
%            = 0 (if didn't lose a cap and all is fine)
%
% NOTES: 
% if crossed a wall, zip it with some probability (p  = 1 if thetacrwall = pi/2)
%  - if it's a long mt --> create a new segment
%  - if it's the beginning of a segment, don't add a new segment, just change the direction. 
%
% ASSUMPTION -- PHYSICS: 
%  The program assumes that the MT cannot bend at an angle that <90 degrees.
%  This means that if a MT approaches (along a wall, for example) a corner, 
%  where the two cell walls intersect at an angle <90, it will have to lose cap
%  and start shrinking. 
%
% Note: here I try to encorporate zipping MTs -- Apr 28, 2015
% Note: implemented program for probzipDM01 based on Deinum and Mulder paper. 
%
% The latest corrections (aug 25) are in the program addonedimer03
% What is new and implemented is marked with (+), the suggestions are marked with (-)
% + zip along the closest one
% - probability of zipping increases if need to pass more MTs
% - marking the crossed ones such that braiding does not happen
% - if the very FIRST dimer of the MT has to zip
%   - make sure that it does not increase the number of segments in a MT 
%    (this might cause the bug that the MT has 2 segments if it starts parallel to the wall,
%   - change the first angle of the MT to be 1/2 in the direction of the MT that is it zipping along. 
%     This way all the MTs will be INSIDE the domain and none of them will grow along the cell-wall from the start. 
%
% >>>>>>>>> the version 4 <<<< <<<<<<<
%  tries to resolve the problem of braiding. Braiding might be biologically
%  similar to what's happening, but in the simulations it is undesirable. 
%  The reason is that it uses up many segments. 
%  Here there are several solutoins to braiding, which are explained in the
%  appropriate sections, e.g. collapse if need to zip twice in one go. 
%


function [site, segm, xseed, intsctysno, intsctdist, intsctdistside,capyesno,bigstate,state] ...
         = proga_addonedimer06(site,siten,segm,sidenode,seedside,xseed,gamma,intsctysno,intsctdist,intsctdistside,N,Nsegments,Nactin,L,intsctinfty,capyesno,bigstate,state,thetacrwall,thetacrmt,Pcat)
     
      
     eps = 1.0e-8; 
     boo_index = 0; 
     
      
     % local variables to this program 
     array_actually_crossed = zeros(N*Nsegments); 
     
     
     segmn = segm(siten);                                  % we're in a particular segment 
     

     site(siten,segmn,1) = site(siten,segmn,1)  + 1;       % grow the MT by 1 in that segment. 
     gamma1 = [0,cumsum(gamma(2:end))];                    % an array of gammas where the angles are computed the same way as the MT angles
     nextsidearr = circshift([1:length(L)]',-1);           % two useful arrays -- the index of the next side and the previous side. 
     prevsidearr = circshift([1:length(L)]',+1); 
     
     % if it's the first dimer, initialize all the xing arrays
     if and(segmn == 1,site(siten,segmn,1) ==1), 
         [intsctdistside,intsctysno,intsctdist] = ...
             proga_init_xing01(N,Nsegments,Nactin,site,siten,segm,seedside,sidenode,xseed,L,gamma,intsctdistside,intsctysno,intsctdist,intsctinfty); 
     end
     
     
     
     
     % ---- if crossing a wall --> zip along it -----
     
     [baddist,badsiden] = min(intsctdistside((siten-1)*Nsegments + segmn,:) - site(siten,segmn,1));  % checking if there are intersections with any cell walls
%==========================================================================================
     if (baddist<=0),                                            % if there are intersections with the first wall
%==========================================================================================
           
            
          %==========================================================================================
          if and(segmn == 1, site(siten,segmn,1) == 1),         % if it is the very first tubulin dimer of a MT, align it along a neighboring wall. 
          %==========================================================================================

              [togrowsegm, thetasegm] = proga_zipyesno01(site(siten,segmn,4),gamma1(badsiden), thetacrwall,100);            % check prob growing along a wall, get the direction along which to grow. 
              if togrowsegm < 0.5,  % if not growing, lose cap, shrink 
                    capyesno(siten) = 0;   bigstate([2*siten-1, 2*siten]) = state(:,3);  % lost cap. :<, back to B0 state. 
                    [site,segm,capyesno,intsctysno,intsctdist,intsctdistside,bigstate] =  ...
                        proga_shrinking001(site,siten,segm,Nsegments,capyesno,intsctysno,intsctdist,intsctdistside,bigstate,state,intsctinfty);
              else % if growig from zero into another direction... 
                   if badsiden == prevsidearr(seedside(siten)),   thetasegm = mod(gamma1(prevsidearr(seedside(siten))) + pi,2*pi);  % if intersecting a side, it's either the next one or the previous one. Pick the angle of the new direction (thetasegm) accordingly. 
                   else                                           thetasegm =     gamma1(nextsidearr(seedside(siten)));
                   end 
                   site(siten,segmn,2:4) = [cos(thetasegm), sin(thetasegm), thetasegm]; 
                   [intsctdistside,intsctysno,intsctdist] =   ...
                       proga_init_xing01(N,Nsegments,Nactin,site,siten,segm,seedside,sidenode,xseed,L,gamma,intsctdistside,intsctysno,intsctdist,intsctinfty); 
                   % >>> XING MT --> IF INTERSECTED MTS ON THE WAY --> CROSS THEM, AND SET THE INTERSECTION TO INTSCTINFTY IN THE XING MTX
                   %  >>> WARNING <<<  --- >>> FOR SOME REASON THE CODE NEVER GETS TO THIS PART.  <<<< 
                   a = (intsctysno(segmn,:)>0.5); locs = locones(a) ; % find if there are intersections and where they happen. locs = an array of indices of the mts that intersect our segmn
                   if length(locs) > 0.5,  % if there were intersections, act as if they have never happened. update array accordingly. 
                                           % scenario -- if go || a wall and cross a bunch of MT's, then shrink back to zero and wait: 
                        capyesno(siten) = 0;   bigstate([2*siten-1, 2*siten]) = state(:,3);  % lost cap. :< 
                        [site,segm,capyesno,intsctysno,intsctdist,intsctdistside,bigstate] = ...
                           proga_shrinking001(site,siten,segm,Nsegments,capyesno,intsctysno,intsctdist,intsctdistside,bigstate,state,intsctinfty);
                   end
              end
         
          
          %==========================================================================================
          else               % if there are intersections with the first wall, but "else" = "it's not the very first dimer" --> create a new segment which goes along the wall
          %==========================================================================================
              site(siten,segmn,1) = site(siten,segmn,1) - 1;  % deleted the segment that we've just grown.
              [togrowsegm, thetasegm] = proga_zipyesno01(site(siten,segmn,4),gamma1(badsiden), thetacrwall,100);
              if togrowsegm < 0.5,                            % if the probability is to lose a cap, do so
                    capyesno(siten) = 0;   bigstate([2*siten-1, 2*siten]) = state(:,1);  % lost cap. :< and shrink (below)
                    [site,segm,capyesno,intsctysno,intsctdist,intsctdistside,bigstate] = ...
                        proga_shrinking001(site,siten,segm,Nsegments,capyesno,intsctysno,intsctdist,intsctdistside,bigstate,state,intsctinfty);
              else                                           %   ['create a new segment on the mt #',num2str(siten)]
                  if segm(siten) + 1 <= Nsegments, 
                      segm(siten) = segm(siten) + 1; % started a new segment
                  else
                      ['GROWING MORE SEGMENTS THAN THERE ARE ALLOWED, (siten,segmn)= ',num2str(siten),', ', num2str(segmn)]
                      pause 
                  end
                  segmn = segm(siten);           % now segmn = the number of the new segment on that MT

                  xseed(siten,segmn,1) = xseed(siten,segmn-1,1) + site(siten,segmn-1,1)*site(siten,segmn-1,2); % create the new seed of the new segment 
                  xseed(siten,segmn,2) = xseed(siten,segmn-1,2) + site(siten,segmn-1,1)*site(siten,segmn-1,3);

                  site(siten,segmn,1)=1; site(siten,segmn,2)=cos(thetasegm); site(siten,segmn,3)=sin(thetasegm); site(siten,segmn,4)=thetasegm;  % update site var's

                  [intsctdistside,intsctysno,intsctdist] = ...
                      proga_init_xing01(N,Nsegments,Nactin,site,siten,segm,seedside,sidenode,xseed,L,gamma,intsctdistside,intsctysno,intsctdist,intsctinfty); 
                  
                  [baddist,badsiden] = min(intsctdistside((siten-1)*Nsegments + segmn,:) - site(siten,segmn,1));   % check if you're not growing straight outside of the domain
                   if (baddist<0),  
                       capyesno(siten) = 0;   bigstate([2*siten-1, 2*siten]) = state(:,1);  % lost cap; 
                       % delete the new segmnt we just started, cause it's no good. 
                       [site,segm,capyesno,intsctysno,intsctdist,intsctdistside,bigstate] = ...
                           proga_shrinking001(site,siten,segm,Nsegments,capyesno,intsctysno,intsctdist,intsctdistside,bigstate,state,intsctinfty);
                       % start shrinking properly. 
                       [site,segm,capyesno,intsctysno,intsctdist,intsctdistside,bigstate] = ...
                           proga_shrinking001(site,siten,segm,Nsegments,capyesno,intsctysno,intsctdist,intsctdistside,bigstate,state,intsctinfty);
                   else 
                       % >>> XING MT --> IF INTERSECTED MTS ON THE WAY --> CROSS THEM, AND SET THE INTERSECTION TO INTSCTINFTY IN THE XING MTX
                       
                       a = (intsctysno(segmn,:)>0.5); locs = locones(a) ; % find if there are intersections and where they happen. locs = an array of indices of the mts that intersect our segmn
                       if length(locs) > 0.5,  %if there were intersections, act as if they have never happened. update array accordingly. 
                        %   ['there are ',num2str(length(locs) ),' intersections with poor MTs - when creating a new segment but not from dimer 1']
                           intsctysno(segmn, locs) = 1;            intsctysno(locs,segmn) = 1;  
                           intsctdist(segmn, locs) = intsctinfty;  intsctdist(locs,segmn) = intsctinfty;                    
                       end  
                       
                       
                       
                   end
              end
          %==========================================================================================
          end   % end the possibility of first crossing the wall -- now should check if crossing a MT
          %==========================================================================================
%==========================================================================================
     else  % not crossing a wall but can cross a MT (!) % might still cross a MT 
%==========================================================================================
         % if not crossing  -> already everything updated, so skip that case
         %
         % If a MT is crossed, then there is a probability of catastrophe, zipping and crossing. 
         % If many MTs are crossed, then if by this algorithm our MT would zip along AT LEAST ONE 
         %    of the MTs, then it zips along the closest one. 
         % If ( the new MT just starts, decide to zip) 
         % AND (theta(MT) along which we zip is the same as gamma1(wall that MT starts from)
         % THEN adjust the angle of the MT
         % 
         
         
         osegmn        = (siten-1)*Nsegments + segmn;  % overall segment number of the segment we're working with
         dist_to_xing  = intsctdist(osegmn,:) - site(siten,segmn,1); 
         locs_prelim   = locones(- dist_to_xing);      % all the MTsegments that we can cross. have to figure out only those who are already there
         locs1         = zeros(N*Nsegments,1); 
         index_crossed = 0;                            % index_crossed = will be the number of MTs that are crossed. 
         siten_crossed_arr = zeros(1,10); 
         segmn_crossed_arr = zeros(1,10);
         togrowsegm_arr    = zeros(1,10);
         thetasegm_arr     = zeros(1,10);
         
         if length(locs_prelim) > 0.5, 
              for i1 = 1:length(locs_prelim), 
                osegmn1 = locs_prelim(i1); 
                
                segmn_crossed         = mod(osegmn1,Nsegments);            % the segment of the MT, that we come to the location of xing to. 
                siten_crossed         = (osegmn1 - segmn_crossed)/Nsegments + 1;
                if segmn_crossed == 0, segmn_crossed = Nsegments; siten_crossed  = siten_crossed  - 1; end      
                
                if intsctdist(osegmn1,osegmn) - site(siten_crossed, segmn_crossed, 1) < 0, 
                    index_crossed = index_crossed + 1; 
                    
                    locs1(index_crossed) = osegmn1; 
                    siten_crossed_arr(index_crossed) = siten_crossed; % remember which sites and segments have crossed
                    segmn_crossed_arr(index_crossed) = segmn_crossed;
                    
                end
             end
         end
                    
         
         if index_crossed >0.5,                                       % if there are interesections
            locs = locs1(1:index_crossed);
             
            [vval,postn] = min(-intsctdist(osegmn,locs));             % find the MT segment with which our MT intersects first. 
            osegmn2 = locs(postn);                                    % overall sgmnt postn of the MT that has the shortest xing with our MT segmnt 
            
            for ii = 1:index_crossed, 
                [togrowsegm_arr(ii), thetasegm_arr(ii)] = ... 
                    proga_zipyesnoDM01(site(siten,segmn,4),site(siten_crossed_arr(ii),segmn_crossed_arr(ii),4),thetacrmt,Pcat);  % check if should zip 
            end 
            thetasegm = thetasegm_arr(postn);                         % in case of zipping, will go along this theta
          
            
            if     prod(togrowsegm_arr (1:index_crossed) - 1) == 0, togrowsegm = 1; % if at least one wants to zip -> zip
            elseif prod(togrowsegm_arr (1:index_crossed)    ) == 0, togrowsegm = 0; % elseif at least one wants to catastrophe -> do so
            else                                                    togrowsegm = 2; % else, cross. 
            end 
            
            
            %*******************
            if  togrowsegm == 0,   %  catastrophe -> togrowsegm = 0, thetasegm = 0; 
                         capyesno(siten) = 0;   bigstate([2*siten-1, 2*siten]) = state(:,1);  % lost cap; delete the new segmnt that just started, no good. 
                         [site,segm,capyesno,intsctysno,intsctdist,intsctdistside,bigstate] = ...
                             proga_shrinking001(site,siten,segm,Nsegments,capyesno,intsctysno,intsctdist,intsctdistside,bigstate,state,intsctinfty);
            %*******************
            elseif togrowsegm == 2,   %  crossing    -> togrowsegm = 2, thetasegm = theta1 --> only update the intersection mtx --> as if nothing happened. 
                intsctysno(osegmn, locs) = 1;            intsctysno(locs,osegmn) = 1;
                intsctdist(osegmn, locs) = intsctinfty;  intsctdist(locs,osegmn) = intsctinfty;  
            %*******************
            else  %if zipping -- zip along the closest one 
                                    
                site(siten,segmn,1) = site(siten,segmn,1) - 1;  % deleted the dimer that we've just grown.
                
                %________________( if it is ***NOT*** the very first dimer, do the following )________________________________
                if not(and( site(siten,segmn,1) == 0 , segm(siten) == 1)), 
                    segm(siten) = segm(siten) + 1;                  % started a new segment only if it's not the very first segment of a MT
                    segmn       = segm(siten);  
                    
                    xseed(siten,segmn,1) = xseed(siten,segmn-1,1) + site(siten,segmn-1,1)*site(siten,segmn-1,2);
                    xseed(siten,segmn,2) = xseed(siten,segmn-1,2) + site(siten,segmn-1,1)*site(siten,segmn-1,3);

                    site(siten,segmn,1)=1; site(siten,segmn,2)=cos(thetasegm); site(siten,segmn,3)=sin(thetasegm); site(siten,segmn,4)=thetasegm;  % update site var's
                    osegmn_new = (siten-1)*Nsegments + segmn;

                    [intsctdistside,intsctysno,intsctdist] = ...
                        proga_init_xing01(N,Nsegments,Nactin,site,siten,segm,seedside,sidenode,xseed,L,gamma,intsctdistside,intsctysno,intsctdist,intsctinfty); 

                    % >> check if it intersects anything now << 
                    
                    [baddist,badsiden] = min(intsctdistside(osegmn_new,:) - site(siten,segmn,1));  % check intersections with walls
                    if (baddist<=eps),    % if intersect a wall --> collapse twice 
                        capyesno(siten) = 0;   bigstate([2*siten-1, 2*siten]) = state(:,1);  % lost cap; delete the just-created dimer
                        [site,segm,capyesno,intsctysno,intsctdist,intsctdistside,bigstate] = ...
                            proga_shrinking001(site,siten,segm,Nsegments,capyesno,intsctysno,intsctdist,intsctdistside,bigstate,state,intsctinfty);
                        % start shrinking properly. 
                        [site,segm,capyesno,intsctysno,intsctdist,intsctdistside,bigstate] = ...
                            proga_shrinking001(site,siten,segm,Nsegments,capyesno,intsctysno,intsctdist,intsctdistside,bigstate,state,intsctinfty);
                    else
                        
                     % check if have crossed any other MTs now that the segment is new. 
                     % if have intersected them, cross them or collapse. 
                     
                     
                             osegmnz       = (siten-1)*Nsegments + segmn;  % overall segment number of the segment we're working with
                             dist_to_xingz  = intsctdist(osegmnz,:) - site(siten,segmn,1); 
                             locs_prelimz   = locones(- dist_to_xingz);      % all the MTsegments that we can cross. have to figure out only those who are already there
                             locs1z         = zeros(N*Nsegments,1); 
                             index_crossedz = 0;                            % index_crossed = will be the number of MTs that are crossed. 
                             siten_crossed_arrz = zeros(1,10); 
                             segmn_crossed_arrz = zeros(1,10);
                             togrowsegm_arrz    = zeros(1,10);
                             thetasegm_arrz     = zeros(1,10);

                             if length(locs_prelimz) > 0.5, 
                                  for i1 = 1:length(locs_prelimz), 
                                    osegmn1z = locs_prelimz(i1); 

                                    segmn_crossedz         = mod(osegmn1z,Nsegments);            % the segment of the MT, that we come to the location of xing to. 
                                    siten_crossedz         = (osegmn1z - segmn_crossedz)/Nsegments + 1;
                                    if segmn_crossedz == 0, segmn_crossedz = Nsegments; siten_crossedz  = siten_crossedz  - 1; end      

                                    if intsctdist(osegmn1z,osegmnz) - site(siten_crossedz, segmn_crossedz, 1) < 0, 
                                        index_crossedz = index_crossedz + 1; 

                                        locs1z(index_crossedz) = osegmn1z; 
                                        siten_crossed_arrz(index_crossedz) = siten_crossedz; % remember which sites and segments have crossed
                                        segmn_crossed_arrz(index_crossedz) = segmn_crossedz;

                                    end
                                 end
                             end


                             if index_crossedz >0.5,                                       % if there are interesections
                                locsz = locs1z(1:index_crossedz);

                                [vval,postnz] = min(-intsctdist(osegmnz,locsz));             % find the MT segment with which our MT intersects first. 
                                osegmn2 = locsz(postnz);                                    % overall sgmnt postn of the MT that has the shortest xing with our MT segmnt 

                                for ii = 1:index_crossedz, 
                                    [togrowsegm_arrz(ii), thetasegm_arrz(ii)] = ... 
                                        proga_zipyesnoDM01(site(siten,segmn,4),site(siten_crossed_arrz(ii),segmn_crossed_arrz(ii),4),thetacrmt,Pcat);  % check if should zip 
                                end 
                                thetasegm = thetasegm_arrz(postnz);                         % in case of zipping, will go along this theta


                                if     prod(togrowsegm_arrz (1:index_crossedz) - 1) == 0, togrowsegm = 2; % if at least one wants to zip -> CROSS 
                                elseif prod(togrowsegm_arrz (1:index_crossedz)    ) == 0, togrowsegm = 0; % elseif at least one wants to catastrophe -> do so
                                else                                                      togrowsegm = 2; % else, cross. 
                                end 


                                %*******************
                                if  togrowsegm == 0,   %  catastrophe -> togrowsegm = 0, thetasegm = 0; 
                                             capyesno(siten) = 0;   bigstate([2*siten-1, 2*siten]) = state(:,1);  % lost cap; delete the new segmnt that just started, no good. 
                                             [site,segm,capyesno,intsctysno,intsctdist,intsctdistside,bigstate] = ...
                                                 proga_shrinking001(site,siten,segm,Nsegments,capyesno,intsctysno,intsctdist,intsctdistside,bigstate,state,intsctinfty);
                                             % and start shrinking properly
                                              [site,segm,capyesno,intsctysno,intsctdist,intsctdistside,bigstate] = ...
                                                 proga_shrinking001(site,siten,segm,Nsegments,capyesno,intsctysno,intsctdist,intsctdistside,bigstate,state,intsctinfty);
                                             
                                %*******************
                                elseif togrowsegm == 2,   %  crossing    -> togrowsegm = 2, thetasegm = theta1 --> only update the intersection mtx --> as if nothing happened. 
                                    intsctysno(osegmnz, locsz) = 1;            intsctysno(locsz,osegmnz) = 1;
                                    intsctdist(osegmnz, locsz) = intsctinfty;  intsctdist(locsz,osegmnz) = intsctinfty;    

                                end 
                             end
                                  
                    end
                %________________( if it ***IS***  the very first dimer... )________________________________
                else  
                    
                    alpha111    = gamma1(seedside(siten)); 
                    segmn       = 1; 
                    segm(siten) = 1; 
                    
                    if or (abs(alpha111 - thetasegm)<1.e-2, abs(alpha111 + pi - thetasegm) < 1.e-2),
                         
                        % figure out if the new growth is parallel to the wall from zero length. 
                        % thetasegm = the direction of the growth of the new segment 
                        % seedside(siten) = the side of the wall that the MT belongs to
                        % gamma1(seedside(siten))
                        site(siten,segmn,1) = 1; 
                        
                        xx0 = xseed(siten,segmn,1);  % seed of the zipping mt
                        yy0 = xseed(siten,segmn,2); 
                        
                        xx1 = xseed(siten_crossed_arr(postn), segmn_crossed_arr(postn),1); % seed of the segment which we zip along
                        yy1 = xseed(siten_crossed_arr(postn), segmn_crossed_arr(postn),2);
                        
                        aa = [xx0 - xx1,yy0 - yy1];
                        bb = [cos(thetasegm),sin(thetasegm)];
                        pp = aa - bb* (aa*bb');
                        dd = bb * sqrt(1 - norm(pp)^2/4) - pp/2; 
                        
                        site(siten,segmn,2) = dd(1); 
                        site(siten,segmn,3) = dd(2);  
                        if dd(2)<0, site(siten,segmn,4) = 2*pi - acos(dd(1));
                        else        site(siten,segmn,4) = acos(dd(1));
                        end
                    
                    else % if the growth is NOT parallel, then do the normal thing
                        % >> DO THIS ONLY IF THE MT DOES NOT START GROWING OUTSIDE THE DOMAIN 
                       if and(thetasegm > gamma1(seedside(siten)), thetasegm < gamma1(seedside(siten)) + pi), 
                            site(siten,segmn,1) = 1; 
                            site(siten,segmn,2) = cos(thetasegm); 
                            site(siten,segmn,3) = sin(thetasegm); 
                            site(siten,segmn,4) = thetasegm;  % update site var's
                       else % if the MT starts zipping outside of the domain, then it has to lose cap and shrink to zero length. 
                           capyesno(siten) = 0;   bigstate([2*siten-1, 2*siten]) = state(:,1);  % lost cap. :< Delete the new dimer, which is the FIRST one.
                           [site,segm,capyesno,intsctysno,intsctdist,intsctdistside,bigstate] = ...
                              proga_shrinking001(site,siten,segm,Nsegments,capyesno,intsctysno,intsctdist,intsctdistside,bigstate,state,intsctinfty);
                           boo_index = 1; 
                          
                          %>>> start growing parallel to the other MT in
                          %the direction of the cell !!!
                          
%                             site(siten,segmn,1) = 1; 
%                             site(siten,segmn,2) = cos(thetasegm+pi); 
%                             site(siten,segmn,3) = sin(thetasegm+pi); 
%                             site(siten,segmn,4) = mod(thetasegm+pi,2*pi); 
%                           
%                           
                       end  
                    end
                    
                    if boo_index < 0.5,  % in case started growing something (i.e. didn't lose cap and shrunk)
                        osegmn_new = (siten-1)*Nsegments + segmn;
                        % in either case -- parallel or NOT , check the intersections with stuff after have just created the new segment. 

                        [intsctdistside,intsctysno,intsctdist] = ...
                            proga_init_xing01(N,Nsegments,Nactin,site,siten,segm,seedside,sidenode,xseed,L,gamma,intsctdistside,intsctysno,intsctdist,intsctinfty); 

                        % >> check if it intersects anything now a. wall, b. MT << 
                        [baddist,badsiden] = min(intsctdistside(osegmn_new,:) - site(siten,segmn,1));  % check intersections with walls
                        if (baddist<=eps),                                                       % if intersect a wall --> collapse   
                            capyesno(siten) = 0;   bigstate([2*siten-1, 2*siten]) = state(:,1);  % lost cap. :< Delete the new doimer, which is the FIRST one.
                            [site,segm,capyesno,intsctysno,intsctdist,intsctdistside,bigstate] = ...
                                proga_shrinking001(site,siten,segm,Nsegments,capyesno,intsctysno,intsctdist,intsctdistside,bigstate,state,intsctinfty);
                        else
                           % there are no intersections with the wall of the new MT segment
                           % if intersect any number of MTs --> pretend like nothing happened
                           % >>> what to put here ??? <<< 
                        end
                    end 
                end
            end
         else % if have reached any of the intersections for the other MTs and the other MTs are not yet there, put that we've been there. 
             
              osegmn = (siten-1)*Nsegments + segmn;                              % overall segment number of the segment we're working with
              dist_to_xing = intsctdist(osegmn,:) - site(siten,segmn,1); 
              locs = locones(-dist_to_xing); 
              if length(locs) > 0.5, 
                  locs = locones(- dist_to_xing .* (1-intsctysno(:,osegmn))');
                  intsctysno(osegmn,locs) = 1; 
              end
         end
      end %__ END OF ZIPPING ALONG MT'S OPTION ____
end
% __ GRAND END OF THE FUNCTION ____ 
  
 
    
    
    