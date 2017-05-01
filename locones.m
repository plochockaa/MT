%========================================================================================== 
% This program finds the locations of all the nonzero elements of a vector
%
% The name of the funciton comes from loc (locations) of ones (ones), therefore locones(a).
% 
% The array (a) DOES NOT HAVE TO BE 0's and 1's, at the first stage it is made to be so. 
% The function can be used to find the positions of the positive elements of an array. 
% 
% INPUT VARIABLES: 
%   a - array of zeros and ones 
% 
% OUTPUT VARIABLES: 
%   locs = array of the positions of positive elements of (a). 
% 
% TO CALL THE FUNCTION: 
%   locs = locones(a); 
%
%==========================================================================================



function locs = locones(a)

    a = a>0; 
    m = sum(a>0.5); 
    
    if m < 0.5,  locs = [];  % in case there are no 1's, locs = []
    elseif m < 1.5, [v, locs] = max(a); % if there is only one 1, locs = position of that 1
    else    % if there are more than one 1, then locs = an array of all the positions of nonzero elements of a
        locs = zeros(1,m);
        [v,locs(1)] = max(a); a(locs(1)) = 0; 
        for i = 2:m, 
            [v,locs(i)] = max(a); 
            a(locs(i)) = 0; 
        end
    end
end