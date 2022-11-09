% it's original V_bus 5*5 matrix
Y_bus = [3.729-j*49.72, 0, 0, 0, -3.729+j*49.72;
         0, 2.679-j*28.46, 0, -0.893+j*9.920, -1.786+j*19.84;
         0, 0, 7.458-j*99.44, -7.458+j*99.44, 0;
         0, -0.893+j*9.920, -7.458+j*99.44, 11.922-j*147.96, -3.571+j*39.68;
         -3.729+j*49.72, -1.786+j*19.84, 0, -3.571+j*39.68, 8.193-j*98.66];

 V = [1; 1 ; 1.05; 1; 1];  % 초기 V값(모선별 V값을 열벡터로 나타내었습니다)
 P = [0; -5; 4.4; 0; 0];  % 초기 P값(미지값(1)은 0으로 나타내었습니다)
 Q = [0; -1.5; 0; 0; 0];  % 초기 Q값(미지값(1)은 0으로 나타내었습니다)
 
 dP = [0; 0; 0; 0];
 dQ = [0; 0; 0; 0];
 
 X_log = [0;0;0;0;1;1.05;1;1]; 
 
 for i = 1:30
     % 1st step) calculate delta(P)
     [P1, Q1] = findPQ(Y_bus, V);
     dP = P(2:5) - P1(2:5);
     dQ = Q(2:5) - Q1(2:5);    
     
     % 2st step) calculate zaccobian matrix
     J = Jacc(Y_bus, V, P, Q);
     
     % 3rd step) calculate delta(X)
     dX = J\[dP;dQ];
     
     % 4th step) calculate next_X
     
     mag = abs(V(2:5))+dX(5:8);
     ang = angle(V(2:5))+dX(1:4);
     V(2:5) = mag.*exp(j * ang);
     X_log = [X_log, X_log(:,i)+dX];
     % 5th step) check converge and calculate error           
     error_count = 0;
     error = [0;0;0;0;0];
     for k = 2:5
         error(k) = abs((X_log(k, i+1) - X_log(k, i))/X_log(k, i));
         if (error(k) <= 10^(-4))
             error_count = error_count + 1;
         end
     end
     % 모두 수렴조건에 만족한다면 for문을 탈출한다.
     if (error_count == 4) 
         break;
     end
 end
 
 
% find real power and reactive power
function [P, Q] = findPQ(Y_bus, V)
     I = [0;0;0;0;0];
    
     I = Y_bus*V;  % current of I
     S = V.*conj(I);  % S = VI*
     P = real(S);     % real factor of S
     Q = imag(S);     % Imag factor of S
end
 
 % find Jaccobian matrix (8*8)
 function [J] = Jacc(Y_bus, V, P, Q)
     J1 = zeros(4, 4);
     J2 = zeros(4, 4);
     J3 = zeros(4, 4);
     J4 = zeros(4, 4);
     for i = 2:5
        for j = 2:5
            if (i == j) % diagonal jacobian element
                for n = 1:5
                    arc = angle(V(i)) - angle(V(n)) - angle(Y_bus(i, n));
                    if (i ~= n)
                        J1(i-1, j-1) = J1(i-1, j-1) -abs(V(i))*abs(Y_bus(i, n))*abs(V(n))*sin(arc);
                        J3(i-1, j-1) = J3(i-1, j-1) +abs(V(i))*abs(Y_bus(i, n))*abs(V(n))*cos(arc);
                    end             
                    J2(i-1, j-1) = J2(i-1, j-1) + abs(Y_bus(i, n))*abs(V(n))*cos(arc);
                    J4(i-1, j-1) = J4(i-1, j-1) + abs(Y_bus(i, n))*abs(V(n))*sin(arc);
                end
                J2(i-1, j-1) = J2(i-1, j-1) + abs(V(i))*abs(Y_bus(i, i))*cos(angle(Y_bus(i, i)));
                J4(i-1, j-1) = J4(i-1, j-1) - abs(V(i))*abs(Y_bus(i, i))*sin(angle(Y_bus(i, i)));
            else  % non-diagonal jacobian element
                arc = angle(V(i)) - angle(V(j)) - angle(Y_bus(i, j));
                J1(i-1, j-1) = abs(V(i))*abs(Y_bus(i, j))*abs(V(j))*sin(arc);
                J2(i-1, j-1) = abs(V(i))*abs(Y_bus(i, j))*cos(arc);
                J3(i-1, j-1) = -abs(V(i))*abs(Y_bus(i, j))*abs(V(j))*cos(arc);
                J4(i-1, j-1) = abs(V(i))*abs(Y_bus(i, j))*sin(arc);
            end
        end
     end
     J = [J1, J2; J3, J4];
 end
