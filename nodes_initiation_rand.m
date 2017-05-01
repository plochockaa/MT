
%==========================================================================
% Description of the nodes_initiation function: 
%==========================================================================
% _____( input )________
% N     = wanted number of seeds
% L     = array(Nsides,1) = lengths of the sides of the cell. The last one does not have to be specified. 
% gamma = array(Nsides,1) = angles between the sides (see drawing), the firs and the last oen are computed in this subroutine. 
%
%  _____( output )________
% xseed    = array(N,Nsegments,2) of N seed sites (x,y) pairs in 2D
%          CORRECTION: the array xseed became the seed array for the seeds
%          of the MT and each of its segments as well (Apr 26, 2015). 
% seedside = array(N,1) = array of which cell side the seed belongs to
% L        = array with all the cell side lengths
% gamma    = array of all the angles between the cell sides
% sidenode = array(Nsides,2) of (x,y) coordinates in 2D of all the cell corners
%
% ___( noes on physical properties )_____-
% - the minimum distance between the nodes is 3, since the tubulin dimer
% height  = 8.2nm  and the outer diameter of the MT is 24nm \approx 3
% heights. 
% - non-dimensionalization for the unit of length is 1 height of the
% tubulin dimer = 8.2nm. 
% 
%==========================================================================



function [xseed,seedside,L,gamma,sidenode,site] = nodes_initiation_rand(N,Nactin,Nsegments,L,gamma)

    Nsides   = length(L); 
    
    sidenode = zeros(Nsides,2);    % nodes = corners of the cell 
   
    sidenode(1,1:2) = [0,   0];    % the first node is at the origin in 2D
    sidenode(2,1:2) = [L(1),0];    % the second node is at [L(1),0] in 2D
    for k = 3:Nsides, 
        sidenode(k,1:2) = sidenode(k-1,1:2) + L(k-1)*[cos(sum(gamma(2:k-1))), sin(sum(gamma(2:k-1)))];
    end
    L(Nsides) = sqrt(sidenode(Nsides,1)^2 + sidenode(Nsides,2)^2); 
    argh = atan(min(sidenode(Nsides,2)/sidenode(Nsides,1),1.e+16)); 
    if argh < 0, argh = argh + pi; end
    gamma(1)      = pi - argh; 
    gamma(Nsides) = 2*pi - sum(gamma(1:Nsides - 1)); 

    %_____( seeding seeding seeding seeding )_____________

    dx = sum(L)/N; 

    xseed    = zeros(N,Nsegments,2);  % (x and y of the seed)
    seedside = zeros(N,1);  % index of the side 1:N that the seed belongs to
    xalong   = dx/2;        % index of how far the seed is along the perimeter
                            % the first point is at dx/2 not to be in the corner
                            % betwen L(1) and L(Nsides). 
    xseed(1,1,1:2) = [xalong,0];  % deal with seed 1 at the origin
    seedside(1)  = 1;    

    for k = 2:N,         %for all the points starting the 2nd one
        xalong = xalong + dx; 

        if xalong < L(1),   % if the seed is on the L1 side  
            seedside(k) = 1; 
            xseed(k,1,1:2) = [xalong,0];
        else
           for m = 1:Nsides - 1, 
               if and( xalong > sum(L(1:m)), xalong < sum(L(1:m+1)) ),
                   seedside(k) = m+1; 
                   xseed(k,1,1) = L(1)  + sum(L(2:m).*cos(cumsum(gamma(2:m)))) + (xalong - sum(L(1:m)))*cos(sum(gamma(2:m+1))); 
                   xseed(k,1,2) =         sum(L(2:m).*sin(cumsum(gamma(2:m)))) + (xalong - sum(L(1:m)))*sin(sum(gamma(2:m+1))); 
               end
           end
        end
    end
    
      
     %%%%%%%%%%% ADD RANDOM ACTIN FILAMENTS %%%%%%%%%%%%
     
     gamma1 = [0,cumsum(gamma(2:end))];
     site     = zeros(N+Nactin,Nsegments,4);  
     
    xv = [sidenode(1,1);sidenode(2,1);sidenode(3,1);sidenode(4,1);sidenode(1,1)];
    yv = [sidenode(1,2);sidenode(2,2);sidenode(3,2);sidenode(4,2);sidenode(1,2)];

    min_x = max(xv);
    max_x = max(xv);
    min_y = max(yv);
    max_y = max(yv);

    pos_actin = zeros(Nactin,2);

    while sum(pos_actin(:,1)==0)>0;

        % generate random points in rectangle with boundaries defined by outer
        % points of trapezium

        if min_x < 0;
            L = max_x + min_x;
            pos_x = rand*L + min_x;
        else
            pos_x = rand*max_x;
        end

        pos_y = rand*max_y;

        [in,on] = inpolygon(pos_x,pos_y,xv,yv);

        if in == 1;

            index = find(pos_actin(:,1)==0);
            pos_actin(index(1),:) = [pos_x pos_y];

        end

    end

    for jj = 1:Nactin;

            r = zeros(1,4);
            x1(jj) = pos_actin(jj,1);
            y1(jj) = pos_actin(jj,2);
            theta1(jj) = rand*2*pi;

            for kk = 1:4
                x2 = sidenode(kk,1);
                y2 = sidenode(kk,2);
                theta2 = gamma1(kk);
                [r1,r2] = intersectionrs01(x1(jj),y1(jj),theta1(jj),x2,y2,theta2);
                r(kk) = (r1<0)*600 + (r1>=0).*r1;
            end

           radius(jj) = min(r);
           x(jj) = radius(jj)*cos(theta1(jj));
           y(jj) = radius(jj)*sin(theta1(jj));



    end

    seeds_actin = [[pos_actin(:,1)+x'], [pos_actin(:,2)+y']];

    for jj = 1:Nactin;


            r = zeros(1,4);
            x1(jj) = seeds_actin(jj,1);
            y1(jj) = seeds_actin(jj,2);
            theta(jj) = mod(theta1(jj)+pi,2*pi);

            for kk = 1:4
                x2 = sidenode(kk,1);
                y2 = sidenode(kk,2);
                theta2 = gamma1(kk);
                [r1,r2] = intersectionrs01(x1(jj),y1(jj),theta(jj),x2,y2,theta2);
                r_new(kk) = (r1<=1.e-5)*600 + (r1>1.e-5).*r1;
            end

           radius_new(jj) = min(r_new);
           x_new(jj) = radius_new(jj)*cos(theta(jj));
           y_new(jj) = radius_new(jj)*sin(theta(jj));

           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    end

    xseed(N+1:N+Nactin,1,1) = seeds_actin(:,1);
    xseed(N+1:N+Nactin,1,2) = seeds_actin(:,2);
    site(N+1:N+Nactin,1,4) = theta;
    site(N+1:N+Nactin,1,2) = seeds_actin(:,1);
    site(N+1:N+Nactin,1,3) = seeds_actin(:,2);
    site(N+1:N+Nactin,1,1) = radius_new;
     
     
     
     
     
 end


