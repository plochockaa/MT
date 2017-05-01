
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



function [actual_Nactin,xseed,seedside,L,gamma,sidenode,site] = nodes_initiation_parallel(N,Nactin,Nsegments,L,gamma,theta_actin,a,ecc_cell)

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

    xseed    = zeros(N+Nactin,Nsegments,2);  % (x and y of the seed)
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
    
    gamma1 = [0,cumsum(gamma(2:end))];
      
     %%%%%%%%%%% ADD PARALLEL ACTIN FILAMENTS %%%%%%%%%%%%
    ecc_standard = 0.98;
    [b_ind] = lengths_for_ecc(ecc_standard,N,gamma,a,Nsegments);
    L_ind = [a b_ind  a b_ind]; 
    [xseed1,seedside1,L1,gamma1,sidenode1] = nodes_initiation(N,Nsegments,L_ind,gamma);

    length1 = sqrt((sidenode1(1,1)-sidenode1(3,1))^2 + (sidenode1(1,2)-sidenode1(3,2))^2);
    length2 = sqrt((sidenode1(2,1)-sidenode1(4,1))^2 + (sidenode1(2,2)-sidenode1(4,2))^2);
    Lcell = max([length1,length2,L1]);

    inc = Lcell/(Nactin+1);
    
    %%%%%%%% NOW WORK WITH THE CELL OF A GIVEN ECC %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [b] = lengths_for_ecc(ecc_cell,N,gamma,a,Nsegments);
    L     = [a b  a b]; 
    [xseed,seedside,L,gamma,sidenode] = nodes_initiation(N,Nsegments,L,gamma);
    gamma1 = [0,cumsum(gamma(2:end))];

    X = [xseed(:,1,1)-mean(xseed(:,1,1)),xseed(:,1,2)-mean(xseed(:,1,2))]; [u,s,v] = svd(X);
    mid_cell = [mean(xseed(:,1,1)),mean(xseed(:,1,2))];

    theta_cell = atan(v(2,1)/v(1,1));

    theta_ACT = theta_cell-(pi/2-theta_actin);
    theta_ACT = mod(theta_ACT,pi);

    % describe the line at the main axis 
    mid_x = mid_cell(1) + 2*Lcell*cos(theta_cell+pi);
    mid_y = mid_cell(2) + 2*Lcell*sin(theta_cell+pi);

    mid_r = 4*Lcell;

    %%%%%%%%%%% FIND POSITIONS OF SIDENODES XING MAIN AXIS %%%%%%%%%%%%%%%%%%%%

    x1 = mid_x;
    y1 = mid_y;
    theta1 = theta_cell;

        for kk = 1:4;

            x2 = sidenode(kk,1);
            y2 = sidenode(kk,2);
            theta2 = theta_ACT;
            [r1,r2] = intersectionrs01(x1,y1,theta1,x2,y2,theta2);
            R(kk) = (r1<0)*0 + (r1>=0).*r1;

        end

        min_R = min(R); % see figure 1 for the idea behind main axis calcs
        max_R = max(R);

        % actual number of actin filaments based on current cell ecc and
        % interval from ecc = 0.98
        actual_Nactin = floor((max_R - min_R)/inc);
        site     = zeros(N+actual_Nactin,Nsegments,4);  
        seed_actin = zeros(actual_Nactin,2);

        for jj = 1:actual_Nactin;

            seed_actin(jj,:) = [mid_x+(min_R+(jj)*inc)*cos(theta_cell),mid_y+(min_R+(jj)*inc)*sin(theta_cell)];

        end

        %%%%%%%%%%%%%%% FIND SEEDS ON SIDE EDGES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % move all the cell positions outside of the cell to find seeds easier!

        seed_actin2 = [seed_actin(:,1)+500*cos(theta_ACT+pi) seed_actin(:,2)+500*sin(theta_ACT+pi)];

        if theta_ACT<pi/2;

            for jj = 1:actual_Nactin;

            x1(jj) = seed_actin2(jj,1);
            y1(jj) = seed_actin2(jj,2);
            theta1 = theta_ACT;

               for kk = 2:3;

                x2 = sidenode(kk,1);
                y2 = sidenode(kk,2);
                theta2 = gamma1(kk);
                [r1,r2] = intersectionrs01(x1(jj),y1(jj),theta1,x2,y2,theta2);
                r(kk) = (r1<1.e-4)*600 + (r1>=1.e-4).*r1;

               end

            site(N+jj,1,2) = x1(jj) + min(r([2 3]))*cos(theta_ACT);
            site(N+jj,1,3) = y1(jj) + min(r([2 3]))*sin(theta_ACT);
            site(N+jj,1,4) = theta_ACT + pi;

            end

        elseif theta_ACT>pi/2;

              for jj = 1:actual_Nactin;

            x1(jj) = seed_actin2(jj,1);
            y1(jj) = seed_actin2(jj,2);
            theta1 = theta_ACT;

               for kk = 3:4;

                x2 = sidenode(kk,1);
                y2 = sidenode(kk,2);
                theta2 = gamma1(kk);
                [r1,r2] = intersectionrs01(x1(jj),y1(jj),theta1,x2,y2,theta2);
                r(kk) = (r1<1.e-4)*600 + (r1>=1.e-4).*r1;

               end

            site(N+jj,1,2) = x1(jj) + min(r([4 3]))*cos(theta_ACT);
            site(N+jj,1,3) = y1(jj) + min(r([4 3]))*sin(theta_ACT);
            site(N+jj,1,4) = theta_ACT + pi;

            end


        end

        xseed(N+1:N+actual_Nactin,1,1) = site(N+1:N+actual_Nactin,1,2);
        xseed(N+1:N+actual_Nactin,1,2) = site(N+1:N+actual_Nactin,1,3);

    %%%%%%%%%%%% actin filaments (FIND RADIUS, we have pos and theta %%%%%%%%%%


      for jj = 1:actual_Nactin;

            x1a = site(N+jj,1,2);
            y1a = site(N+jj,1,3);
            theta1 = site(N+jj,1,4);

             for kk = 1:4;

                x2 = sidenode(kk,1);
                y2 = sidenode(kk,2);
                theta2 = gamma1(kk);
                [r1,r2] = intersectionrs01(x1a,y1a,theta1,x2,y2,theta2);
                r(kk) = (r1<1.e-4)*600 + (r1>=1.e-4)*r1;

             end

            if min(r)<600;

                site(N+jj,1,1) = min(r);

            else

                site(N+jj,1,1) = 0;

            end


      end
     
     
     
     
 end


