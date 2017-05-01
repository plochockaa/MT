% shrinking code: 
% 
% 1. have to shrink the segment
% 2. if the segment disappears, segm(siten) = segm(siten) - 1; 
%    --> update the intersection matrix 
% 3. if the MT shrinks to zero, it gains a cap
% 4. if the MT shrinks to zero lengths, it has ZERO segments. 
% 
% ________________________________________________________________________
% To call the shinking code from the main program: 
% 
%  [site,segm,capyesno,intsctysno,intsctdist,intsctdistside,bigstate] = ...
%         proga_shrinking001 ...
%         (site,segm,Nsegments,capyesno,intsctysno,intsctdist,intsctdistside,bigstate,state)
% ________________________________________________________________________

function [site,segm,capyesno,intsctysno,intsctdist,intsctdistside,bigstate] = ...
         proga_shrinking001(site,siten,segm,Nsegments,capyesno,intsctysno,intsctdist,intsctdistside,bigstate,state,intsctinfty)

tubn = (siten-1)*Nsegments + segm(siten);               % intro the number of the tubulin dimers in the MT

site(siten,segm(siten),1) = site(siten,segm(siten),1)  - 1;    % shrink


if  site(siten,segm(siten),1) < 0.5,                 % if length(segment) = 0, it is not a segment
	% update the intersection matrix                 % update the intersection matrices
    intsctysno(tubn,:)     = 0; 
    intsctysno(:,tubn)     = 0; 
    intsctdist(tubn,:)     = intsctinfty;
    intsctdist(:,tubn)     = intsctinfty;
    intsctdistside(tubn,:) = intsctinfty;
    
    segm(siten) = segm(siten)-1;                     % delete the segment 
end 
if site(siten,1,1) <0.5,                             %      if shrunk to zero MT length
    bigstate([2*siten-1, 2*siten]) = state(:,3);     %          change the bigstate to be able to only gain caps
    intsctysno(tubn,:)     = 0;                     %          update the intersectino matrices.
    intsctysno(:,tubn)     = 0;                     %
    intsctdist(tubn,:)     = intsctinfty;           %
    intsctdist(:,tubn)     = intsctinfty;           %
    intsctdistside(tubn,:) = intsctinfty;           %                                                                                     
end                                                  %      end


 


