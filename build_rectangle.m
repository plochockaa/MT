%_________________________________________________________________
% build_rectangle.m builds a rectangle around the given MT segment. 
% 
% INPUT: 
%    seed_corrds(1,2) = column vector of the two coordinates of the seed of the segment
%    r          = length of the segment
%    cos_theta  = cos(theta) of the segment
%    sin_theta  = sin(theta) of the segment 
%    half_width = half-width of the rectangle = 1.5 for the MTs for the 
%
% OUTPUT
%    A = [p1, p2, p3, p4], where  p's are 2D column vector coord's 
%                        of the corners of the rectangle starting with the
%                        "bottom left" one from the seed and going in the 
%                        counterclockwise direction.  
%_________________________________________________________________

function A = build_rectangle(seed_coords, r, cos_theta, sin_theta, half_width)

    n = [-sin_theta, cos_theta]'; % normal vector to the MT segment 
    l = [ cos_theta, sin_theta]'; % direction vector of the MT segment 
    
    A = zeros(2,4); 
    
    A(:,1) = seed_coords - n*half_width;
    A(:,2) = A(:,1) + r*l; 
    A(:,4) = A(:,1) + 2*n*half_width; 
    A(:,3) = A(:,4) + r*l; 
    
end