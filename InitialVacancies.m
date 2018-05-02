clear;

%We will simulate the idea for treating stripe formation as a result of an
%ASEP-like process in a 2D grid with two strains.  Each site will have four values:  2,1,
%-1,-2. 1 corresponds to a particle site being occupied by a vertically
%oriented cell; -1 is analogously defined for horizontally oriented cells.
%A site having |1| means it is of strain 1.  2 corresponds to a particle
%site being occupied by a vertically oriented cell; -2 is defined
%analogously for horizontally oriented cells.  A site having |2| means it
%is occupied by a cell of the second strain. 0 indicates a particular site is unoccupied.  We will run a Monte-Carlo
%simulation and look at the ending result.  

%Parameter initialization
Lx  = .1;
Ly  = 1;
eps = .01;
N   = Lx/eps;
M   = Ly/eps;
h   =  5;
a   =  10;
b   =  10;
dt  = .01;
C   = 5000000;
aa  = ones(4,1);
lambda = 0;
kappa = 1;

%Initialization of grid
numbers = [2,1,0,-1,-2];

A = zeros(M,N);
for i = 1:N
    for j = 1:M
        number = randi([1,length(numbers)]);
        A(j,i) = numbers(number);
    end
end
% A(:,1:N/2) = 1;
% A(:,N/2+1:end) = 2;
%A(M/2,(N/2-1:N/2+1)) = -1;

%Rules
%Cell will grow and divide only in the direction in which it is oriented.
%Cells do not move on their own.  If site adjacent to cell in the direction of its preferred orientation is vacant,
%growth will cause vacant site to become occupied by a cell of the same
%orientation.  If they are occupied, then growth will occur but cells will
%push the ones in front or behind by one unit of the fundamental length
%scale.

for j = 1:C
    R1 = randi([1,M]);
    R2 = randi([1,N]);
    sum1 = 0;
    sum2 = 0;
    sum3 = 0;
    check = A(R1,R2);
    if(R1>1 && R1<M && R2>1 && R2<N) 
        if(check == 1 || check == 2)
            up1 = R1;
            up2 = R1;
            while(up1+1 <= M && abs(A(up1+1,R2))>0)
                sum1 = sum1 + 1;
                up1  = up1 + 1;
            end
            while(up2-1 >= 1 && abs(A(up2-1,R2))>0)
                sum2 = sum2 + 1;
                up2  = up2 - 1;
            end
            for k = R1
                for l = [R2-1, R2+1]
                    if(sign(A(k,l)) == sign(check))
                        sum3 = sum3 + 1;
                    end
                end
            end
            for k = [R1-1,R1+1]
                for l = R2
                    if(sign(A(k,l)) == sign(check))
                        sum3 = sum3 + 1;
                    end
                end
            end
            aa(1) = h/(1+kappa*sum1)*dt;
            aa(2) = h/(1+kappa*sum2)*dt;
            aa(3) = (a*(4-sum3)+lambda)*dt;
            aa(4) = 1 - aa(1) - aa(2)-aa(3);
           
            selection_intv1 = cumsum(aa); 
            selection_intv1 = selection_intv1./selection_intv1(end);
    
            ind = find(selection_intv1 > rand,1);
            
            if(ind == 1 && up1 == R1)
                A(R1 + 1, R2) = A(R1,R2);
            elseif(ind == 2 && up2 == R1)
                A(R1 - 1, R2) = A(R1, R2);
            elseif(ind == 1 && up1 == M)
                for k = up1:-1:R1+1
                    A(k,R2) = A(k-1,R2);
                end
            elseif(ind ==2 && up2 == 1)
                for k = up2:1:R1-1
                    A(k,R2) = A(k+1,R2);
                end
            elseif(ind == 1)
                for k = up1+1:-1:R1+1
                    A(k,R2) = A(k-1,R2);
                end
            elseif(ind == 2)
                for k = up2-1:1:R1-1
                    A(k,R2) = A(k+1,R2);
                end
            elseif(ind == 3)
                A(R1,R2) = -A(R1,R2);
            end
            
        elseif(check == -1 || check == -2)
            up1 = R2;
            up2 = R2;
            while(up1 + 1 <= N && abs(A(R1,up1+1))>0)
                sum1 = sum1 + 1;
                up1  = up1 + 1;
            end
            while(up2 - 1 >= 1 && abs(A(R1,up2-1))>0)
                sum2 = sum2 + 1;
                up2  = up2 - 1;
            end
            for k = R1
                for l = [R2-1, R2+1]
                    if(sign(A(k,l)) == sign(check))
                        sum3 = sum3 + 1;
                    end
                end
            end
            for k = [R1-1,R1+1]
                for l = R2
                    if(sign(A(k,l)) == sign(check))
                        sum3 = sum3 + 1;
                    end
                end
            end
            aa(1) = h/(1+kappa*sum1)*dt;
            aa(2) = h/(1+kappa*sum2)*dt;
            aa(3) = (b*(4-sum3)+lambda)*dt;
            aa(4) = 1-aa(1)-aa(2)-aa(3);
    
            selection_intv1 = cumsum(aa); 
            selection_intv1 = selection_intv1./selection_intv1(end);
    
            ind = find(selection_intv1 > rand,1);
            
            if(ind == 1 && up1 == R2)
                A(R1, R2+1) = A(R1,R2);
            elseif(ind == 2 && up2 == R2)
                A(R1, R2-1) = A(R1,R2);
            elseif(ind == 1 && up1 == N)
                for k = up1:-1:R2+1
                    A(R1,k) = A(R1,k-1);
                end
            elseif(ind ==2 && up2 == 1)
                for k = up2:1:R2-1
                    A(R1,k) = A(R1,k+1);
                end
            elseif(ind == 1)
                for k = up1+1:-1:R2+1
                    A(R1,k) = A(R1, k-1);
                end
            elseif(ind == 2)
                for k = up2-1:1:R2-1
                    A(R1,k) = A(R1, k+1);
                end
            elseif(ind == 3)
                A(R1,R2) = -A(R1,R2);
            end
        end
    elseif(R1 == M && R2 ~= 1 && R2 ~= N)
        if(check == 1 || check == 2)
            up1 = R1;
            up2 = R1;
%             while(up1+1 <= M && abs(A(up1+1,R2))>0)
%                 sum1 = sum1 + 1;
%                 up1  = up1 + 1;
%             end
            while(up2-1 >= 1 && abs(A(up2-1,R2))>0)
                sum2 = sum2 + 1;
                up2  = up2 - 1;
            end
            for k = R1
                for l = [R2-1, R2+1]
                    if(sign(A(k,l)) == sign(check))
                        sum3 = sum3 + 1;
                    end
                end
            end
            for k = R1-1
                for l = R2
                    if(sign(A(k,l)) == sign(check))
                        sum3 = sum3 + 1;
                    end
                end
            end
            aa(1) = 0;%h/(1+kappa*sum1)*dt;
            aa(2) = h/(1+kappa*sum2)*dt;
            aa(3) = (a*(3-sum3)+lambda)*dt;
            aa(4) = 1-aa(1)-aa(2)-aa(3);
           
            selection_intv1 = cumsum(aa); 
            selection_intv1 = selection_intv1./selection_intv1(end);
    
            ind = find(selection_intv1 > rand,1);
            
            if(ind == 1 && up1 == R1)
                A(R1 + 1, R2) = A(R1,R2);
            elseif(ind == 2 && up2 == R1)
                A(R1 - 1, R2) = A(R1, R2);
            elseif(ind == 1 && up1 == M)
                for k = up1:-1:R1+1
                    A(k,R2) = A(k-1,R2);
                end
            elseif(ind ==2 && up2 == 1)
                for k = up2:1:R1-1
                    A(k,R2) = A(k+1,R2);
                end
            elseif(ind == 1)
                for k = up1+1:-1:R1+1
                    A(k,R2) = A(k-1,R2);
                end
            elseif(ind == 2)
                for k = up2-1:1:R1-1
                    A(k,R2) = A(k+1,R2);
                end
            elseif(ind == 3)
                A(R1,R2) = -A(R1,R2);
            end
            
        elseif(check == -1 || check == -2)
            up1 = R2;
            up2 = R2;
            while(up1 + 1 <= N && abs(A(R1,up1+1))>0)
                sum1 = sum1 + 1;
                up1  = up1 + 1;
            end
            while(up2 - 1 >= 1 && abs(A(R1,up2-1))>0)
                sum2 = sum2 + 1;
                up2  = up2 - 1;
            end
            for k = R1
                for l = [R2-1, R2+1]
                    if(sign(A(k,l)) == sign(check))
                        sum3 = sum3 + 1;
                    end
                end
            end
            for k = R1-1
                for l = R2
                    if(sign(A(k,l)) == sign(check))
                        sum3 = sum3 + 1;
                    end
                end
            end
            aa(1) = h/(1+kappa*sum1)*dt;
            aa(2) = h/(1+kappa*sum2)*dt;
            aa(3) = (b*(3-sum3)+lambda)*dt;
            aa(4) = 1-aa(1)-aa(2)-aa(3);
    
            selection_intv1 = cumsum(aa); 
            selection_intv1 = selection_intv1./selection_intv1(end);
    
            ind = find(selection_intv1 > rand,1);
            
            if(ind == 1 && up1 == R2)
                A(R1, R2+1) = A(R1,R2);
            elseif(ind == 2 && up2 == R2)
                A(R1, R2-1) = A(R1,R2);
            elseif(ind == 1 && up1 == N)
                for k = up1:-1:R2+1
                    A(R1,k) = A(R1,k-1);
                end
            elseif(ind ==2 && up2 == 1)
                for k = up2:1:R2-1
                    A(R1,k) = A(R1,k+1);
                end
            elseif(ind == 1)
                for k = up1+1:-1:R2+1
                    A(R1,k) = A(R1, k-1);
                end
            elseif(ind == 2)
                for k = up2-1:1:R2-1
                    A(R1,k) = A(R1, k+1);
                end
            elseif(ind == 3)
                A(R1,R2) = -A(R1,R2);
            end
        end
    elseif(R1 == 1 && R2 ~= 1 && R2 ~= N)
        if(check == 1 || check == 2)
            up1 = R1;
            up2 = R1;
            while(up1+1 <= M && abs(A(up1+1,R2))>0)
                sum1 = sum1 + 1;
                up1  = up1 + 1;
            end
%             while(up2-1 >= 1 && abs(A(up2-1,R2))>0)
%                 sum2 = sum2 + 1;
%                 up2  = up2 - 1;
%             end
            for k = R1
                for l = [R2-1, R2+1]
                    if(sign(A(k,l)) == sign(check))
                        sum3 = sum3 + 1;
                    end
                end
            end
            for k = R1+1
                for l = R2
                    if(sign(A(k,l)) == sign(check))
                        sum3 = sum3 + 1;
                    end
                end
            end
            aa(1) = h/(1+kappa*sum1)*dt;
            aa(2) = 0;%h/(1+kappa*sum2)*dt;
            aa(3) = (a*(3-sum3)+lambda)*dt;
            aa(4) = 1-aa(1)-aa(2)-aa(3);
           
            selection_intv1 = cumsum(aa); 
            selection_intv1 = selection_intv1./selection_intv1(end);
    
            ind = find(selection_intv1 > rand,1);
            
            if(ind == 1 && up1 == R1)
                A(R1 + 1, R2) = A(R1,R2);
            elseif(ind == 2 && up2 == R1)
                A(R1 - 1, R2) = A(R1, R2);
            elseif(ind == 1 && up1 == M)
                for k = up1:-1:R1+1
                    A(k,R2) = A(k-1,R2);
                end
            elseif(ind ==2 && up2 == 1)
                for k = up2:1:R1-1
                    A(k,R2) = A(k+1,R2);
                end
            elseif(ind == 1)
                for k = up1+1:-1:R1+1
                    A(k,R2) = A(k-1,R2);
                end
            elseif(ind == 2)
                for k = up2-1:1:R1-1
                    A(k,R2) = A(k+1,R2);
                end
            elseif(ind == 3)
                A(R1,R2) = -A(R1,R2);
            end
            
        elseif(check == -1 || check == -2)
            up1 = R2;
            up2 = R2;
            while(up1 + 1 <= N && abs(A(R1,up1+1))>0)
                sum1 = sum1 + 1;
                up1  = up1 + 1;
            end
            while(up2 - 1 >= 1 && abs(A(R1,up2-1))>0)
                sum2 = sum2 + 1;
                up2  = up2 - 1;
            end
            for k = R1
                for l = [R2-1, R2+1]
                    if(sign(A(k,l)) == sign(check))
                        sum3 = sum3 + 1;
                    end
                end
            end
            for k = R1+1
                for l = R2
                    if(sign(A(k,l)) == sign(check))
                        sum3 = sum3 + 1;
                    end
                end
            end
            aa(1) = h/(1+kappa*sum1)*dt;
            aa(2) = h/(1+kappa*sum2)*dt;
            aa(3) = (b*(3-sum3)+lambda)*dt;
            aa(4) = 1-aa(1)-aa(2)-aa(3);
    
            selection_intv1 = cumsum(aa); 
            selection_intv1 = selection_intv1./selection_intv1(end);
    
            ind = find(selection_intv1 > rand,1);
            
            if(ind == 1 && up1 == R2)
                A(R1, R2+1) = A(R1,R2);
            elseif(ind == 2 && up2 == R2)
                A(R1, R2-1) = A(R1,R2);
            elseif(ind == 1 && up1 == N)
                for k = up1:-1:R2+1
                    A(R1,k) = A(R1,k-1);
                end
            elseif(ind ==2 && up2 == 1)
                for k = up2:1:R2-1
                    A(R1,k) = A(R1,k+1);
                end
            elseif(ind == 1)
                for k = up1+1:-1:R2+1
                    A(R1,k) = A(R1, k-1);
                end
            elseif(ind == 2)
                for k = up2-1:1:R2-1
                    A(R1,k) = A(R1, k+1);
                end
            elseif(ind == 3)
                A(R1,R2) = -A(R1,R2);
            end
        end
    elseif(R2 == N && R1 ~= 1 && R1 ~= M)
        if(check == 1 || check == 2)
            up1 = R1;
            up2 = R1;
            while(up1+1 <= M && abs(A(up1+1,R2))>0)
                sum1 = sum1 + 1;
                up1  = up1 + 1;
            end
            while(up2-1 >= 1 && abs(A(up2-1,R2))>0)
                sum2 = sum2 + 1;
                up2  = up2 - 1;
            end
            for k = R1
                for l = R2-1
                    if(sign(A(k,l)) == sign(check))
                        sum3 = sum3 + 1;
                    end
                end
            end
            for k = [R1-1,R1+1]
                for l = R2
                    if(sign(A(k,l)) == sign(check))
                        sum3 = sum3 + 1;
                    end
                end
            end
            aa(1) = h/(1+kappa*sum1)*dt;
            aa(2) = h/(1+kappa*sum2)*dt;
            aa(3) = (a*(3-sum3)+lambda)*dt;
            aa(4) = 1-aa(1)-aa(2)-aa(3);
           
            selection_intv1 = cumsum(aa); 
            selection_intv1 = selection_intv1./selection_intv1(end);
    
            ind = find(selection_intv1 > rand,1);
            
            if(ind == 1 && up1 == R1)
                A(R1 + 1, R2) = A(R1,R2);
            elseif(ind == 2 && up2 == R1)
                A(R1 - 1, R2) = A(R1, R2);
            elseif(ind == 1 && up1 == M)
                for k = up1:-1:R1+1
                    A(k,R2) = A(k-1,R2);
                end
            elseif(ind ==2 && up2 == 1)
                for k = up2:1:R1-1
                    A(k,R2) = A(k+1,R2);
                end
            elseif(ind == 1)
                for k = up1+1:-1:R1+1
                    A(k,R2) = A(k-1,R2);
                end
            elseif(ind == 2)
                for k = up2-1:1:R1-1
                    A(k,R2) = A(k+1,R2);
                end
            elseif(ind == 3)
                A(R1,R2) = -A(R1,R2);
            end
            
        elseif(check == -1 || check == -2)
            up1 = R2;
            up2 = R2;
%             while(up1 + 1 <= N && abs(A(R1,up1+1))>0)
%                 sum1 = sum1 + 1;
%                 up1  = up1 + 1;
%             end
            while(up2 - 1 >= 1 && abs(A(R1,up2-1))>0)
                sum2 = sum2 + 1;
                up2  = up2 - 1;
            end
            for k = R1
                for l = R2-1
                    if(sign(A(k,l)) == sign(check))
                        sum3 = sum3 + 1;
                    end
                end
            end
            for k = [R1-1,R1+1]
                for l = R2
                    if(sign(A(k,l)) == sign(check))
                        sum3 = sum3 + 1;
                    end
                end
            end
            aa(1) = 0;%h/(1+kappa*sum1)*dt;
            aa(2) = h/(1+kappa*sum2)*dt;
            aa(3) = (b*(3-sum3)+lambda)*dt;
            aa(4) = 1-aa(1)-aa(2)-aa(3);
    
            selection_intv1 = cumsum(aa); 
            selection_intv1 = selection_intv1./selection_intv1(end);
    
            ind = find(selection_intv1 > rand,1);
            
            if(ind == 1 && up1 == R2)
                A(R1, R2+1) = A(R1,R2);
            elseif(ind == 2 && up2 == R2)
                A(R1, R2-1) = A(R1,R2);
            elseif(ind == 1 && up1 == N)
                for k = up1:-1:R2+1
                    A(R1,k) = A(R1,k-1);
                end
            elseif(ind ==2 && up2 == 1)
                for k = up2:1:R2-1
                    A(R1,k) = A(R1,k+1);
                end
            elseif(ind == 1)
                for k = up1+1:-1:R2+1
                    A(R1,k) = A(R1, k-1);
                end
            elseif(ind == 2)
                for k = up2-1:1:R2-1
                    A(R1,k) = A(R1, k+1);
                end
            elseif(ind == 3)
                A(R1,R2) = -A(R1,R2);
            end
        end
    elseif(R2 == 1 && R1 ~= 1 && R1 ~= M)
        if(check == 1 || check == 2)
            up1 = R1;
            up2 = R1;
            while(up1+1 <= M && abs(A(up1+1,R2))>0)
                sum1 = sum1 + 1;
                up1  = up1 + 1;
            end
            while(up2-1 >= 1 && abs(A(up2-1,R2))>0)
                sum2 = sum2 + 1;
                up2  = up2 - 1;
            end
            for k = R1
                for l = R2+1
                    if(sign(A(k,l)) == sign(check))
                        sum3 = sum3 + 1;
                    end
                end
            end
            for k = [R1-1,R1+1]
                for l = R2
                    if(sign(A(k,l)) == sign(check))
                        sum3 = sum3 + 1;
                    end
                end
            end
            aa(1) = h/(1+kappa*sum1)*dt;
            aa(2) = h/(1+kappa*sum2)*dt;
            aa(3) = (a*(3-sum3)+lambda)*dt;
            aa(4) = 1-aa(1)-aa(2)-aa(3);
           
            selection_intv1 = cumsum(aa); 
            selection_intv1 = selection_intv1./selection_intv1(end);
    
            ind = find(selection_intv1 > rand,1);
            
            if(ind == 1 && up1 == R1)
                A(R1 + 1, R2) = A(R1,R2);
            elseif(ind == 2 && up2 == R1)
                A(R1 - 1, R2) = A(R1, R2);
            elseif(ind == 1 && up1 == M)
                for k = up1:-1:R1+1
                    A(k,R2) = A(k-1,R2);
                end
            elseif(ind ==2 && up2 == 1)
                for k = up2:1:R1-1
                    A(k,R2) = A(k+1,R2);
                end
            elseif(ind == 1)
                for k = up1+1:-1:R1+1
                    A(k,R2) = A(k-1,R2);
                end
            elseif(ind == 2)
                for k = up2-1:1:R1-1
                    A(k,R2) = A(k+1,R2);
                end
            elseif(ind == 3)
                A(R1,R2) = -A(R1,R2);
            end
            
        elseif(check == -1 || check == -2)
            up1 = R2;
            up2 = R2;
            while(up1 + 1 <= N && abs(A(R1,up1+1))>0)
                sum1 = sum1 + 1;
                up1  = up1 + 1;
            end
%             while(up2 - 1 >= 1 && abs(A(R1,up2-1))>0)
%                 sum2 = sum2 + 1;
%                 up2  = up2 - 1;
%             end
            for k = R1
                for l = R2+1
                    if(sign(A(k,l)) == sign(check))
                        sum3 = sum3 + 1;
                    end
                end
            end
            for k = [R1-1,R1+1]
                for l = R2
                    if(sign(A(k,l)) == sign(check))
                        sum3 = sum3 + 1;
                    end
                end
            end
            aa(1) = h/(1+kappa*sum1)*dt;
            aa(2) = 0;%h/(1+kappa*sum2)*dt;
            aa(3) = (b*(3-sum3)+lambda)*dt;
            aa(4) = 1-aa(1)-aa(2)-aa(3);
    
            selection_intv1 = cumsum(aa); 
            selection_intv1 = selection_intv1./selection_intv1(end);
    
            ind = find(selection_intv1 > rand,1);
            
            if(ind == 1 && up1 == R2)
                A(R1, R2+1) = A(R1,R2);
            elseif(ind == 2 && up2 == R2)
                A(R1, R2-1) = A(R1,R2);
            elseif(ind == 1 && up1 == N)
                for k = up1:-1:R2+1
                    A(R1,k) = A(R1,k-1);
                end
            elseif(ind ==2 && up2 == 1)
                for k = up2:1:R2-1
                    A(R1,k) = A(R1,k+1);
                end
            elseif(ind == 1)
                for k = up1+1:-1:R2+1
                    A(R1,k) = A(R1, k-1);
                end
            elseif(ind == 2)
                for k = up2-1:1:R2-1
                    A(R1,k) = A(R1, k+1);
                end
            elseif(ind == 3)
                A(R1,R2) = -A(R1,R2);
            end
        end
    elseif(R1 == 1 && R2 == 1)
        if(check == 1 || check == 2)
            up1 = R1;
            up2 = R1;
            while(up1+1 <= M && abs(A(up1+1,R2))>0)
                sum1 = sum1 + 1;
                up1  = up1 + 1;
            end
%             while(up2-1 >= 1 && abs(A(up2-1,R2))>0)
%                 sum2 = sum2 + 1;
%                 up2  = up2 - 1;
%             end
            for k = R1
                for l = R2+1
                    if(sign(A(k,l)) == sign(check))
                        sum3 = sum3 + 1;
                    end
                end
            end
            for k = R1+1
                for l = R2
                    if(sign(A(k,l)) == sign(check))
                        sum3 = sum3 + 1;
                    end
                end
            end
            aa(1) = h/(1+kappa*sum1)*dt;
            aa(2) = 0;%h/(1+kappa*sum2)*dt;
            aa(3) = (a*(2-sum3)+lambda)*dt;
            aa(4) = 1-aa(1)-aa(2)-aa(3);
           
            selection_intv1 = cumsum(aa); 
            selection_intv1 = selection_intv1./selection_intv1(end);
    
            ind = find(selection_intv1 > rand,1);
            
            if(ind == 1 && up1 == R1)
                A(R1 + 1, R2) = A(R1,R2);
            elseif(ind == 2 && up2 == R1)
                A(R1 - 1, R2) = A(R1, R2);
            elseif(ind == 1 && up1 == M)
                for k = up1:-1:R1+1
                    A(k,R2) = A(k-1,R2);
                end
            elseif(ind ==2 && up2 == 1)
                for k = up2:1:R1-1
                    A(k,R2) = A(k+1,R2);
                end
            elseif(ind == 1)
                for k = up1+1:-1:R1+1
                    A(k,R2) = A(k-1,R2);
                end
            elseif(ind == 2)
                for k = up2-1:1:R1-1
                    A(k,R2) = A(k+1,R2);
                end
            elseif(ind == 3)
                A(R1,R2) = -A(R1,R2);
            end
            
        elseif(check == -1 || check == -2)
            up1 = R2;
            up2 = R2;
            while(up1 + 1 <= N && abs(A(R1,up1+1))>0)
                sum1 = sum1 + 1;
                up1  = up1 + 1;
            end
%             while(up2 - 1 >= 1 && abs(A(R1,up2-1))>0)
%                 sum2 = sum2 + 1;
%                 up2  = up2 - 1;
%             end
            for k = R1
                for l = R2+1
                    if(sign(A(k,l)) == sign(check))
                        sum3 = sum3 + 1;
                    end
                end
            end
            for k = R1+1
                for l = R2
                    if(sign(A(k,l)) == sign(check))
                        sum3 = sum3 + 1;
                    end
                end
            end
            aa(1) = h/(1+kappa*sum1)*dt;
            aa(2) = 0;%h/(1+kappa*sum2)*dt;
            aa(3) = (b*(2-sum3)+lambda)*dt;
            aa(4) = 1-aa(1)-aa(2)-aa(3);
    
            selection_intv1 = cumsum(aa); 
            selection_intv1 = selection_intv1./selection_intv1(end);
    
            ind = find(selection_intv1 > rand,1);
            
            if(ind == 1 && up1 == R2)
                A(R1, R2+1) = A(R1,R2);
            elseif(ind == 2 && up2 == R2)
                A(R1, R2-1) = A(R1,R2);
            elseif(ind == 1 && up1 == N)
                for k = up1:-1:R2+1
                    A(R1,k) = A(R1,k-1);
                end
            elseif(ind ==2 && up2 == 1)
                for k = up2:1:R2-1
                    A(R1,k) = A(R1,k+1);
                end
            elseif(ind == 1)
                for k = up1+1:-1:R2+1
                    A(R1,k) = A(R1, k-1);
                end
            elseif(ind == 2)
                for k = up2-1:1:R2-1
                    A(R1,k) = A(R1, k+1);
                end
            elseif(ind == 3)
                A(R1,R2) = -A(R1,R2);
            end
        end
    elseif(R1 == 1 && R2 == N)
        if(check == 1 || check == 2)
            up1 = R1;
            up2 = R1;
            while(up1+1 <= M && abs(A(up1+1,R2))>0)
                sum1 = sum1 + 1;
                up1  = up1 + 1;
            end
%             while(up2-1 >= 1 && abs(A(up2-1,R2))>0)
%                 sum2 = sum2 + 1;
%                 up2  = up2 - 1;
%             end
            for k = R1
                for l = R2-1
                    if(sign(A(k,l)) == sign(check))
                        sum3 = sum3 + 1;
                    end
                end
            end
            for k = R1+1
                for l = R2
                    if(sign(A(k,l)) == sign(check))
                        sum3 = sum3 + 1;
                    end
                end
            end
            aa(1) = h/(1+kappa*sum1)*dt;
            aa(2) = 0;%h/(1+kappa*sum2)*dt;
            aa(3) = (a*(2-sum3)+lambda)*dt;
            aa(4) = 1-aa(1)-aa(2)-aa(3);
           
            selection_intv1 = cumsum(aa); 
            selection_intv1 = selection_intv1./selection_intv1(end);
    
            ind = find(selection_intv1 > rand,1);
            
            if(ind == 1 && up1 == R1)
                A(R1 + 1, R2) = A(R1,R2);
            elseif(ind == 2 && up2 == R1)
                A(R1 - 1, R2) = A(R1, R2);
            elseif(ind == 1 && up1 == M)
                for k = up1:-1:R1+1
                    A(k,R2) = A(k-1,R2);
                end
            elseif(ind ==2 && up2 == 1)
                for k = up2:1:R2-1
                    A(k,R2) = A(k+1,R2);
                end
            elseif(ind == 1)
                for k = up1+1:-1:R1+1
                    A(k,R2) = A(k-1,R2);
                end
            elseif(ind == 2)
                for k = up2-1:1:R1-1
                    A(k,R2) = A(k+1,R2);
                end
            elseif(ind == 3)
                A(R1,R2) = -A(R1,R2);
            end
            
        elseif(check == -1 || check == -2)
            up1 = R2;
            up2 = R2;
%             while(up1 + 1 <= N && abs(A(R1,up1+1))>0)
%                 sum1 = sum1 + 1;
%                 up1  = up1 + 1;
%             end
            while(up2 - 1 >= 1 && abs(A(R1,up2-1))>0)
                sum2 = sum2 + 1;
                up2  = up2 - 1;
            end
            for k = R1
                for l = R2-1
                    if(sign(A(k,l)) == sign(check))
                        sum3 = sum3 + 1;
                    end
                end
            end
            for k = R1+1
                for l = R2
                    if(sign(A(k,l)) == sign(check))
                        sum3 = sum3 + 1;
                    end
                end
            end
            aa(1) = 0;%h/(1+kappa*sum1)*dt;
            aa(2) = h/(1+kappa*sum2)*dt;
            aa(3) = (b*(2-sum3)+lambda)*dt;
            aa(4) = 1-aa(1)-aa(2)-aa(3);
    
            selection_intv1 = cumsum(aa); 
            selection_intv1 = selection_intv1./selection_intv1(end);
    
            ind = find(selection_intv1 > rand,1);
            
            if(ind == 1 && up1 == R2)
                A(R1, R2+1) = A(R1,R2);
            elseif(ind == 2 && up2 == R2)
                A(R1, R2-1) = A(R1,R2);
            elseif(ind == 1 && up1 == N)
                for k = up1:-1:R1+1
                    A(R1,k) = A(R1,k-1);
                end
            elseif(ind ==2 && up2 == 1)
                for k = up2:1:R2-1
                    A(R1,k) = A(R1,k+1);
                end
            elseif(ind == 1)
                for k = up1+1:-1:R2+1
                    A(R1,k) = A(R1, k-1);
                end
            elseif(ind == 2)
                for k = up2-1:1:R2-1
                    A(R1,k) = A(R1, k+1);
                end
            elseif(ind == 3)
                A(R1,R2) = -A(R1,R2);
            end
        end
    elseif(R1 == M && R2 == 1)
        if(check == 1 || check == 2)
            up1 = R1;
            up2 = R1;
%             while(up1+1 <= M && abs(A(up1+1,R2))>0)
%                 sum1 = sum1 + 1;
%                 up1  = up1 + 1;
%             end
            while(up2-1 >= 1 && abs(A(up2-1,R2))>0)
                sum2 = sum2 + 1;
                up2  = up2 - 1;
            end
            for k = R1
                for l = R2+1
                    if(sign(A(k,l)) == sign(check))
                        sum3 = sum3 + 1;
                    end
                end
            end
            for k = R1-1
                for l = R2
                    if(sign(A(k,l)) == sign(check))
                        sum3 = sum3 + 1;
                    end
                end
            end
            aa(1) = 0;%h/(1+kappa*sum1)*dt;
            aa(2) = h/(1+kappa*sum2)*dt;
            aa(3) = (a*(2-sum3)+lambda)*dt;
            aa(4) = 1-aa(1)-aa(2)-aa(3);
           
            selection_intv1 = cumsum(aa); 
            selection_intv1 = selection_intv1./selection_intv1(end);
    
            ind = find(selection_intv1 > rand,1);
            
            if(ind == 1 && up1 == R1)
                A(R1 + 1, R2) = A(R1,R2);
            elseif(ind == 2 && up2 == R1)
                A(R1 - 1, R2) = A(R1, R2);
            elseif(ind == 1 && up1 == M)
                for k = up1:-1:R1+1
                    A(k,R2) = A(k-1,R2);
                end
            elseif(ind ==2 && up2 == 1)
                for k = up2:1:R1-1
                    A(k,R2) = A(k+1,R2);
                end
            elseif(ind == 1)
                for k = up1+1:-1:R1+1
                    A(k,R2) = A(k-1,R2);
                end
            elseif(ind == 2)
                for k = up2-1:1:R1-1
                    A(k,R2) = A(k+1,R2);
                end
            elseif(ind == 3)
                A(R1,R2) = -A(R1,R2);
            end
            
        elseif(check == -1 || check == -2)
            up1 = R2;
            up2 = R2;
            while(up1 + 1 <= N && abs(A(R1,up1+1))>0)
                sum1 = sum1 + 1;
                up1  = up1 + 1;
            end
%             while(up2 - 1 >= 1 && abs(A(R1,up2-1))>0)
%                 sum2 = sum2 + 1;
%                 up2  = up2 - 1;
%             end
            for k = R1
                for l = R2+1
                    if(sign(A(k,l)) == sign(check))
                        sum3 = sum3 + 1;
                    end
                end
            end
            for k = R1-1
                for l = R2
                    if(sign(A(k,l)) == sign(check))
                        sum3 = sum3 + 1;
                    end
                end
            end
            aa(1) = h/(1+kappa*sum1)*dt;
            aa(2) = 0;%h/(1+kappa*sum2)*dt;
            aa(3) = (b*(2-sum3)+lambda)*dt;
            aa(4) = 1-aa(1)-aa(2)-aa(3);
    
            selection_intv1 = cumsum(aa); 
            selection_intv1 = selection_intv1./selection_intv1(end);
    
            ind = find(selection_intv1 > rand,1);
            
            if(ind == 1 && up1 == R2)
                A(R1, R2+1) = A(R1,R2);
            elseif(ind == 2 && up2 == R2)
                A(R1, R2-1) = A(R1,R2);
            elseif(ind == 1 && up1 == N)
                for k = up1:-1:R2+1
                    A(R1,k) = A(R1,k-1);
                end
            elseif(ind ==2 && up2 == 1)
                for k = up2:1:R2-1
                    A(R1,k) = A(R1,k+1);
                end
            elseif(ind == 1)
                for k = up1+1:-1:R2+1
                    A(R1,k) = A(R1, k-1);
                end
            elseif(ind == 2)
                for k = up2-1:1:R2-1
                    A(R1,k) = A(R1, k+1);
                end
            elseif(ind == 3)
                A(R1,R2) = -A(R1,R2);
            end
        end
    else
        if(check == 1 || check == 2)
            up1 = R1;
            up2 = R1;
%             while(up1+1 <= M && abs(A(up1+1,R2))>0)
%                 sum1 = sum1 + 1;
%                 up1  = up1 + 1;
%             end
            while(up2-1 >= 1 && abs(A(up2-1,R2))>0)
                sum2 = sum2 + 1;
                up2  = up2 - 1;
            end
            for k = R1
                for l = R2-1
                    if(sign(A(k,l)) == sign(check))
                        sum3 = sum3 + 1;
                    end
                end
            end
            for k = R1-1
                for l = R2
                    if(sign(A(k,l)) == sign(check))
                        sum3 = sum3 + 1;
                    end
                end
            end
            aa(1) = 0;%h/(1+kappa*sum1)*dt;
            aa(2) = h/(1+kappa*sum2)*dt;
            aa(3) = (a*(2-sum3)+lambda)*dt;
            aa(4) = 1-aa(1)-aa(2)-aa(3);
           
            selection_intv1 = cumsum(aa); 
            selection_intv1 = selection_intv1./selection_intv1(end);
    
            ind = find(selection_intv1 > rand,1);
            
            if(ind == 1 && up1 == R1)
                A(R1 + 1, R2) = A(R1,R2);
            elseif(ind == 2 && up2 == R1)
                A(R1 - 1, R2) = A(R1, R2);
            elseif(ind == 1 && up1 == M)
                for k = up1:-1:R1+1
                    A(k,R2) = A(k-1,R2);
                end
            elseif(ind ==2 && up2 == 1)
                for k = up2:1:R1-1
                    A(k,R2) = A(k+1,R2);
                end
            elseif(ind == 1)
                for k = up1+1:-1:R1+1
                    A(k,R2) = A(k-1,R2);
                end
            elseif(ind == 2)
                for k = up2-1:1:R1-1
                    A(k,R2) = A(k+1,R2);
                end
            elseif(ind == 3)
                A(R1,R2) = -A(R1,R2);
            end
            
        elseif(check == -1 || check == -2)
            up1 = R2;
            up2 = R2;
%             while(up1 + 1 <= N && abs(A(R1,up1+1))>0)
%                 sum1 = sum1 + 1;
%                 up1  = up1 + 1;
%             end
            while(up2 - 1 >= 1 && abs(A(R1,up2-1))>0)
                sum2 = sum2 + 1;
                up2  = up2 - 1;
            end
            for k = R1
                for l = R2-1
                    if(sign(A(k,l)) == sign(check))
                        sum3 = sum3 + 1;
                    end
                end
            end
            for k = R1-1
                for l = R2
                    if(sign(A(k,l)) == sign(check))
                        sum3 = sum3 + 1;
                    end
                end
            end
            aa(1) = 0;%h/(1+kappa*sum1)*dt;
            aa(2) = h/(1+kappa*sum2)*dt;
            aa(3) = (b*(2-sum3)+lambda)*dt;
            aa(4) = 1-aa(1)-aa(2)-aa(3);
    
            selection_intv1 = cumsum(aa); 
            selection_intv1 = selection_intv1./selection_intv1(end);
    
            ind = find(selection_intv1 > rand,1);
            
            if(ind == 1 && up1 == R2)
                A(R1, R2+1) = A(R1,R2);
            elseif(ind == 2 && up2 == R2)
                A(R1, R2-1) = A(R1,R2);
            elseif(ind == 1 && up1 == N)
                for k = up1:-1:R2+1
                    A(R1,k) = A(R1,k-1);
                end
            elseif(ind ==2 && up2 == 1)
                for k = up2:1:R2-1
                    A(R1,k) = A(R1,k+1);
                end
            elseif(ind == 1)
                for k = up1+1:-1:R2+1
                    A(R1,k) = A(R1, k-1);
                end
            elseif(ind == 2)
                for k = up2-1:1:R2-1
                    A(R1,k) = A(R1, k+1);
                end
            elseif(ind == 3)
                A(R1,R2) = -A(R1,R2);
            end
        end
    end
    if(mod(j,10000)==0)
        figure(1)
        imagesc(A)
        colormap(jet(5))
        colorbar
        set(gca,'fontsize',20)
        drawnow
    end
% %pause(.1)
% %axis([1 100 1 100])
% pause
end
B = zeros(2*M,2*N);
B(1:M,1:N) = A;
cc = -2;
for j = M+1:1:M+5
    B(j,:) = cc;
    cc = cc+1;
end
figure(1)
imagesc(B)
colormap(jet(5))
colorbar
set(gca,'fontsize',20)
xlabel('\mum')
ylabel('\mum')
% z = legend('V2','V1','Vacant','H1','H2');
% set(z,'box','off')
axis([1 N 1 M])

