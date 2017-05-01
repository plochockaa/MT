% function [togrowsegm, thetasegm] = proga_zipyesno01(theta1,theta2,thetacr,paramprob)
%
% finds if in case of potential zipping, one should zip or not. 
% In case of zipping, it affirms that one should grow a segment 
% (togrowsegm = 1) and gives the direction of the new growth (thetasegm). 
% In case of no zipping, it affirms that one should not grow 
% (togrowsegm = 0) and puts a default 0 as the direction of the
% non-existent "new" segment (thetasegm = 0); 
%
% PARAMETERS:
%   theta1 = theta of the MT segment that we are growing
%   theta2 = theta of the MT/wall along which we might zip
%   thetacr = critical value of theta for zipping or not
%   paramprob = parameter for the sigmoid zip/nozip function that gives the
%   thresholds for zipping probabailistically. Normally paramprob = 100. 


function [togrowsegm, thetasegm] = proga_zipyesno01(theta1,theta2, thetacr,paramprob)

   th3 = acos(cos(theta1-theta2));
   th3 = (th3<pi/2)*th3+(th3>pi/2)*(pi-th3); 
   
   if rand < probcross01(th3,thetacr,paramprob), 
       togrowsegm = 1; 
       if cos(theta2-theta1) > 0, thetasegm = theta2; 
       else                       thetasegm = mod(theta2-pi,2*pi); 
       end
   else
       togrowsegm = 0; thetasegm = 0; 
   end
   
end