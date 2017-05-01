function [B] = lengths_for_ecc(ecc_array,N,gamma,a,Nsegments);

B_og = 0*ecc_array;
B    = 0*ecc_array;
t_max  = 10000;
ecc_new2 = zeros(1,t_max);
eps = 1e-3;
eps2 = eps^2;

    for jj = 1:length(ecc_array);
            ECC2 = ecc_array(jj)^2;

            change = zeros(1,1000000);

            t = 1;
            change(1) = 0;
            inc = 0.5;
            b   = a./sqrt(1 - ECC2); 
            B_og(jj) = b;
            L     = [a b  a b]; 


            [xseed,seedside,L,gamma,sidenode] = nodes_initiation(N,Nsegments,L,gamma);

            X = [xseed(:,1,1)-mean(xseed(:,1,1)), xseed(:,1,2)-mean(xseed(:,1,2))];
            [u,s,v] = svd(X);


            ecc_new2(t) =  1 - (s(2,2)/s(1,1))^2;


            while (abs(ecc_new2(t)-ECC2)>eps2) && (t < t_max);


                if ecc_new2(t) > ECC2;
                    t = t+1;
                    change(t) = -1;
                    b = b - inc;

                else
                    t = t+1;
                    change(t) = 1;
                    b = b + inc;
                end 


                L     = [a b a b]; 

                [xseed,seedside,L,gamma,sidenode] = nodes_initiation(N,Nsegments,L,gamma);

                X = [xseed(:,1,1), xseed(:,1,2)];
                X = [xseed(:,1,1)-mean(xseed(:,1,1)), xseed(:,1,2)-mean(xseed(:,1,2))];
                [u,s,v] = svd(X);
                ecc_new2(t) =  1- (s(2,2)/s(1,1))^2;

                if change(t)~=change(t-1);
                    inc = inc/2;
                end



            B(jj) = b;


            end


    end

end
