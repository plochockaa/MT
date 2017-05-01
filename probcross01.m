% FUNCTION probcross01(theta, thetacr, param)
% probcross (version 01) gives the probability of zippering. 
% the funciton in the verison  01 is tanh. 
% 
% PARAMETERS: 
% param gives the choice of how steep tanh is at thetacr (critical value of
% theta) 
% 
% INTENDED USE: 
% if a random draw is smaller than probcross, then zipping happens, 
% if the draw is larger than the function, then MT goes through. 
%
% FUNCTION CHOICE: 
% the function is intended for the beginning to be close to a jump
% function, so tanh is chosen. 
%
% THOUGHTS: 
% the program is implemented so that there is a one-parameter family of
% funcitons that gives probability of zipping/crossing. This choice reduces
% the number of parameters to fit this function to one. 
%
% EXAMPLE OF USE: 
% qwer = probcross01(theta, thetacr, param)
% qwer = probcross01(theta, pi/6,    100)
%
%=================================================================

function qwer = probcross01(theta, thetacr, param)

    t1 = tanh(param*thetacr); 
    t2 = tanh(param*(thetacr - pi/2)); 
    qwer =  1./(t1 - t2).*(- t2 + tanh(param*(thetacr - theta))); 

end