% function qwer = intersectionrs01 (x1,y1,theta1,x2,y2,theta2)
% 
% is the function for finidng the r's to intersections (version 01) from
% two pints (x1,y1) and (x2,y2) in the directions theta1 and theta2. 
% 


function  [r1,r2] = intersectionrs01(x1,y1,theta1,x2,y2,theta2)

    s1  = sin(theta1); s2 = sin(theta2); c1 = cos(theta1); c2 = cos(theta2);
    det = s1*c2 - c1*s2; 
    r1 = (-s2*(x2-x1) + c2*(y2-y1))./det; 
    r2 = (-s1*(x2-x1) + c1*(y2-y1))./det; 
    
end



%_____________________________________________________________________
%     % test of the subroutine: 
%
%     figure(1); clf; hold on; 
%     plot(x1,y1,'or'); plot(x2,y2,'or'); 
%     a = [0:0.01:5]; 
%     plot(x1 + cos(theta1)*a,y1 + sin(theta1)*a,'-b'); 
%     plot(x2 + cos(theta2)*a,y2 + sin(theta2)*a,'-g'); 
%     plot(x1 + cos(theta1)*r1,y1 + sin(theta1)*r1,'db');
%     plot(x2 + cos(theta2)*r2,y2 + sin(theta2)*r2,'om','MarkerSize',10);
%     