% function randarchonside gives (x,y) unitary directions for the MT that grows from a specific cell side numbered (nside). 
% gamma  = array of the angles between the sides. 
% theta1 = the angle with respect to the horizontal. Not redundant, since will need to average it for the real problem. 

function qwer = randarchonside(Nside,gamma)

    gamma1 = [0,cumsum(gamma(2:end))];
   
   % th  = mod(pi*rand  + sum(gamma(2:Nside)),2*pi); 
    th  = mod(pi*rand  + gamma1(Nside),2*pi); 
    x   = cos(th); 
    y   = sin(th); 
    
    qwer = [x,y,th];
end