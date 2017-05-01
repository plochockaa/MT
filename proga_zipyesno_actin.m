% function  [togrowsegm, thetasegm] = proga_zipyesnoDM01(theta1,theta2, thetacr,Pcat)
%
% finds if in case of potential zipping, one should zip or not. 
% the procedure is based on the behavior in plants, paper by 
% Modelling the role of microtubules in cell morphology Deinum and Mulder
%
% PARAMETERS:
%   theta1  = theta of the MT segment that we are growing
%   theta2  = theta of the MT/wall along which we might zip
%   thetacr = critical value of theta for zip versus cross
%   Pcat    = probability of catastrophe
%
% OUTPUT: 
%  catastrophe -> togrowsegm = 0, thetasegm = 0;      -- catastr
%  zipping     -> togrowsegm = 1, thetasegm = th3     -- zip
%  crossing    -> togrowsegm = 2, thetasegm = theta1  -- cross

function [togrowsegm, thetasegm] = proga_zipyesno_actin(pcat_actin,theta1,theta2)

   th3 = acos(cos(theta1-theta2));
   th3 = (th3<pi/2)*th3+(th3>pi/2)*(pi-th3); 
   
     
       if rand < (th3/pi)*pcat_actin,    
           togrowsegm = 0; 
           thetasegm  = 0; % catastrophe
       else
           togrowsegm = 2;                % cross
           thetasegm = theta1; 
           
       end 
   
end
    