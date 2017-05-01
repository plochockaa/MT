%_______________________________________________________________
% this function is given two segments and checks if they are "nearby" or not
%
% INPUT: 
%   ( segment 1 )
%   seed_corrds1 = seed coords of the 1st segment, 2D column vector  
%   seed_coords2 = seed coords of the 2nd segment, 2D column vector 
%   r, cos_theta, sin_theta = for both segments 
%   half_width = half-width of the rectangle. Normally = 1.5
% 
%   (the two functions used are): 
%         qwer = proga_point_inside_or_not(r,A)
%         A    = build_rectangle(seed_coords, r, cos_theta, sin_theta, half_width)
%         qwerqwer = proga_point_inside_or_not_2(y,seed_coords,r,cos_theta,sin_theta,half_width)
% OUTPUT: 
%   index_inside = an indicator if the two recangles have an overlapping section or not.
%                = 0 if no overlap
%                = 1 if there is an overlap
% 
% IMPORTANT NOTE: 
%   the case when the two rectangles intersect at a large angle and no
%   corners of a rectangle is inside the other one is not considered. 
%   The bio explanation of this is that then the MTs intersect, and the
%   brightness does not increase in that case.
%
% RECTANGLES vs QUADRILATERALS
%   The code proga_point_inside_or_not works for any quadrilaterals
%   The code proga_point_inside_or_not_2 works only for rectangles
%_______________________________________________________________


function index_inside = nearby_or_not(seed_coords_1, r1, cos_theta1, sin_theta1, ...
                                      seed_coords_2, r2, cos_theta2, sin_theta2, ...
                                      half_width); 

    index_inside = 0; 
    
    % build both rectangles 
    A1 = build_rectangle(seed_coords_1, r1, cos_theta1, sin_theta1, half_width);
    A2 = build_rectangle(seed_coords_2, r2, cos_theta2, sin_theta2, half_width);
    
    % check wrt to the first rectangle
    for jj = 1:4, 
        if index_inside <0.5,  
          %  index_inside = index_inside + proga_point_inside_or_not(A2(:,jj), A1);
          index_inside = index_inside + proga_point_inside_or_not_2(A2(:,jj), seed_coords_1,r1,cos_theta1,sin_theta1,half_width);
        end
    end
    
    % check wrt second rectangle 
    for jj = 1:4,   
        if index_inside < 0.5, 
           % index_inside = index_inside + proga_point_inside_or_not(A1(:,jj), A2);
            index_inside = index_inside + proga_point_inside_or_not_2(A1(:,jj), seed_coords_2,r2,cos_theta2,sin_theta2,half_width);
        end
    end


end

                                   