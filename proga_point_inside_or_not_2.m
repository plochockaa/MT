%_______________________________________________________________
% this function determines if a point is inside the rectangle or not
%
% INPUT : 
%    y           = 2D column vector = coord's of the point inside/outside
%    r           = 1D lenght of the MT segment 
%    seed_corrds = 2D column vector of the seed of the MT segment 
%    cos_theta, sin_theta = direction of the MT segment 
%    half_width  = half width of the rectangle
%    
% 
% OUTPUT: 
%    1 - the point is inside
%    0 - the point is outside 
%
% NOTE: 
%    this code is faster than its previous versoin
%    proga_point_inside_or_not.m, but it only works for rectangles. 
%    The previous code works for any convex chetyrehugol'nik. : ) 
%
%_______________________________________________________________


function qwerqwer = proga_point_inside_or_not_2(y,seed_coords,r,cos_theta,sin_theta,half_width)

  % aa = coords of the point "inside" the rectangle in the coord system of the rectangle 
  aa = [cos_theta, sin_theta; - sin_theta, cos_theta]*(y - seed_coords - [cos_theta, sin_theta]'*half_width);
  
  if and( and(aa(1)>0,aa(2)>0), and(aa(1)<r,aa(2)<2*half_width)), 
      qwerqwer = 1; 
  else
      qwerqwer = 0; 
  end
  
  
end


