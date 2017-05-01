%========================================================================
% initialize intersection arrays with walls and other MT segments. 
%
% for every segment, have xseed
% for every wall, have sidenode for the seec of the wall. 
%
% INPUT VARIABLES: 
%  intsctysno = (N*Nsegmn, N*Nsegmn)
%  intsctdist = (N*Nsegmn, N*Nsegmn)
%  intsctdistside = (N*Nsegmn,Nsides) - r of the MT segment to the wall. 
% 
% OUTPUT VARIABLES: 
%  intsctysno = (N*Nsegmn, N*Nsegmn)
%  intsctdist = (N*Nsegmn, N*Nsegmn)
%  intsctdistside = (N*Nsegmn,Nsides) - r of the MT segment to the wall. 
%
% CONVENTIONS: 
%  - All the entries in the intsctdist and intsctdistside are non-negative. 
%  - If the distance to xing is too large, or negative --> it's stored as intsctinfty. 
%  - Assume that the array (site) already has the information about the new segment
% 
%========================================================================
 
function [intsctdistside,intsctysno,intsctdist] = ...
         proga_init_xing01(N,Nsegments,Nactin,site,siten,segm,seedside,sidenode,xseed,L,gamma,intsctdistside,intsctysno,intsctdist,intsctinfty)
 
    segmn = segm(siten);
     
    %___ initialize the arrays to something sensible first____
    segmn_overall = (siten-1)*Nsegments + segmn;
    
    intsctdist(segmn_overall,:)  = intsctinfty;   intsctdist(:,segmn_overall)  = intsctinfty; 
    intsctysno(segmn_overall,:)  = 0;             intsctysno(:,segmn_overall)  = 0; 
    intsctdistside(segmn_overall,:) = intsctinfty; 
     
    %______ check interesections with cell walls__________
     
    gamma1 = [0,cumsum(gamma(2:end))];        % angles of the cell sides wrt horizontal  

    if and(site(siten,segmn,1)==1,segm(siten)==1),
          tocheckside = circshift([1:length(L)]',-seedside(siten))';  tocheckside = tocheckside(1:end-1); 
    else  tocheckside = [1:length(L)];
    end
    
    
    for jj = 1:length(tocheckside),      % going through all the cell sides except the one that MT starts at
        walln = tocheckside(jj);         % the number of the wall that we're checking xing with. 
       
        x1 = xseed(siten,segmn,1); y1 = xseed(siten,segmn,2);  theta1  = site(siten,segmn,4);  % index 1 <-- the new MT segment that we're considering. 
        x2 = sidenode(walln,1);    y2 = sidenode(walln,2);     theta2  = gamma1(walln);        % index 2 <-- the WALL that the MT can intersect with. 

        [r,r2] = intersectionrs01(x1,y1,theta1,x2,y2,theta2);
        r = (r<0)*intsctinfty + (r>=0).*r;  r = min(intsctinfty,r);
        intsctdistside(segmn_overall,walln) = r;
    end   %___end( check interesections with cell walls )__________
        
     %______ check interesections with other MT segments __________
     tocheckMTarray = circshift([1:N+Nactin]',-siten)'; tocheckMTarray = tocheckMTarray(1:end-1); % the array of which MTs to check xing with 
     
     for jj = 1:N-1, 
         tocheckMT = tocheckMTarray(jj);
        
         segm(N+1:N+Nactin)=1;
         
         if segm(tocheckMT)>0.5, % if that MT has any segments at all --> check intersections with all the existing segments. 
             for tochecksegm = 1:segm(tocheckMT), % go through all the segments that MT has 
                 
                 % the next line is already implemented in the loop above about xing MT segment and walls. 
                 %[x1,y1] = xseed(siten,segmn,:);  theta1  = site(siten,segmn,4);  % index 1 <-- the new MT segment that we're considering. 
                 x2  = xseed(tocheckMT,tochecksegm,1);  y2 = xseed(tocheckMT,tochecksegm,2); theta2  = site(tocheckMT,tochecksegm,4);  % 2 = index of the segment we're checking
                 [r1,r2] = intersectionrs01(x1,y1,theta1,x2,y2,theta2);
                 tocheck_segmn_overall = (tocheckMT-1)*Nsegments + tochecksegm; % the number of the checked segment overall, not just along its MT
                 
                 if or(r1<0, r2<0), 
                     intsctysno(segmn_overall, tocheck_segmn_overall) = 0;            intsctysno(tocheck_segmn_overall,segmn_overall) = 0; 
                     intsctdist(segmn_overall, tocheck_segmn_overall) = intsctinfty;  intsctdist(tocheck_segmn_overall,segmn_overall) = intsctinfty;
                 else
                     if   r1 <= 1,  intsctysno(segmn_overall, tocheck_segmn_overall) = 1;
                     else          intsctysno(segmn_overall, tocheck_segmn_overall) = 0; 
                     end
                     if   r2 <= 1,  intsctysno(tocheck_segmn_overall,segmn_overall) = 1; 
                     else          intsctysno(tocheck_segmn_overall,segmn_overall) = 0;  
                     end
                     intsctdist(segmn_overall,tocheck_segmn_overall) = min(r1,intsctinfty); 
                     intsctdist(tocheck_segmn_overall,segmn_overall) = min(r2,intsctinfty);
                 end 
             end
         end
         
     end %______ end( check interesections with other MT segments )__________
end %_____( grand end )_____________