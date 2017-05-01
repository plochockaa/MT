function [angle_distributions] = find_angles(transform_seed,transform_xings)
% transform_seed - point considered
% transform_xings 

    for jj = 1:length(transform_xings);
        
        if  transform_seed(1) <= transform_xings(1,jj); % (0,pi/2)

        angle_distributions(jj) = atan((transform_xings(2,jj)-transform_seed(2))/(transform_xings(1,jj)-transform_seed(1)));

        else   % (pi/2,pi)

        angle_distributions(jj) = pi - atan((transform_xings(2,jj)-transform_seed(2))/(transform_seed(1)-transform_xings(1,jj)));

        end
        
    end

end

