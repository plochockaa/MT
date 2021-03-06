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

function [togrowsegm, thetasegm] = proga_zipyesnoDM01(theta1,theta2, thetacr,Pcat)

   th3 = acos(cos(theta1-theta2));
   th3 = (th3<pi/2)*th3+(th3>pi/2)*(pi-th3); 
   
   if th3 < thetacr,            % if theta < thetacr
       if rand < th3*Pcat/thetacr,    
           togrowsegm = 0; 
           thetasegm  = 0; % catastrophe
       else
           togrowsegm = 1;                % zip
           if cos(theta2-theta1) > 0, thetasegm = theta2; 
           else                       thetasegm = mod(theta2-pi,2*pi); 
           end
       end 
   else                            % else (if theta >= thetacr)
       if rand < Pcat,   togrowsegm = 0;  thetasegm = 0;        % if prob < Pcat -> catastrophe 
       else              togrowsegm = 2;  thetasegm = theta1;   % else --> cross
       end
   end
end
    