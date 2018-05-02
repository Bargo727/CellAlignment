clear;

%Here we will run a mean field approximation of the spatial moran model for
%a single strain with two possible orientations.  There will be no
%switching between states, i.e. there is no reduction to an Ising model.

%N is the number of vertical lattice sites and M is the number of
%horizontal lattice sites.
dk = 0.001;
K  = (0:dk:1)';
ll = length(K);
EE = zeros(ll,1);
MM = zeros(ll,1);
FF = zeros(ll,1);
WW = zeros(ll,1);

AR = zeros(ll,1);

N = 30;
M = 10;

for w = 1:ll
    w
kappa = K(w);

%Initialization of the grid
n1  = zeros(M,N);
%n2 = 1-n1;

%Setting up the hopping rates.  We note that by taking $\kappa$ = 0, the
%preferential nature of the hopping will be eliminated
hp = zeros(M,1);
hm = zeros(M,1); 
gp = zeros(N,1);
gm = zeros(N,1);

%Basal hopping rate
h0 = 5;

%Creating the hopping matrices
for i = 1:M
    hp(i) = h0*exp(-kappa*(M-i));
    hm(i) = h0*exp(-kappa*(i-1));
end

for j = 1:N
    gp(j) = h0*exp(-kappa*(N-j));
    gm(j) = h0*exp(-kappa*(j-1));
end

m1 = 1-n1;

F  = zeros(M*N);
B = zeros(N);
BU = zeros(N,1);
BL = zeros(N,1);
for i = 1:M
    K1 = zeros(M);
    K2 = zeros(M);
    K3 = zeros(M);
    K1(i,i) = 1;
    hp = h0*exp(-kappa*(M-i));
    hm = h0*exp(-kappa*(i-1));
    for j = 1:M
        if((j-i)==1)
            K2(i,j) = 1;
        elseif((i-j)==1)
            K3(i,j) = 1;
        end
    end
    
    if(i == 1)
        B(1,1)   = -hm*n1(i+1,1) - gm(2)*(1-n1(i,2));
        B(1,2)   = gm(2)*n1(i,1);
        B(N,N)   = -hm*n1(i+1,N)- gp(N-1)*(1-n1(i,N-1));
        B(N,N-1) = gp(N-1)*n1(i,N);
        
        for j = 2:N-1
            for k = 2:N-1
                if(j==k)
                    B(j,k) = -hm*n1(i+1,j)-gp(k-1)*(1-n1(i,j-1))-gm(k+1)*(1-n1(i,j+1));
                elseif(j-k == 1)
                    B(j,k) = gp(k)*n1(i,j);
                elseif(j-k == -1)
                    B(j,k) = gm(k)*n1(i,j);
                end
            end
        end
    
        for j = 1:N
            BU(j) = hm*(1-n1(i,j));
            BL(j) = hp*(1-n1(i,j));
        end
        
    elseif(i == M)
        B(1,1)   = -hp*n1(i-1,1) - gm(2)*(1-n1(i,2));
        B(1,2)   = gm(2)*n1(i,1);
        B(N,N)   = -hp*n1(i-1,N)- gp(N-1)*(1-n1(i,N-1));
        B(N,N-1) = gp(N-1)*n1(i,N);
        
        for j = 2:N-1
            for k = 2:N-1
                if(j==k)
                    B(j,k) = -hp*n1(i-1,j)-gp(k-1)*(1-n1(i,j-1))-gm(k+1)*(1-n1(i,j+1));
                elseif(j-k == 1)
                    B(j,k) = gp(k)*n1(i,j);
                elseif(j-k == -1)
                    B(j,k) = gm(k)*n1(i,j);
                end
            end
        end
    
        for j = 1:N
            BU(j) = hm*(1-n1(i,j));
            BL(j) = hp*(1-n1(i,j));
        end
        
    else
        B(1,1)   = -hp*n1(i-1,1)-hm*n1(i+1,1) - gm(2)*(1-n1(i,2));
        B(1,2)   = gm(2)*n1(i,1);
        B(N,N)   = -hp*n1(i-1,1)-hm*n1(i+1,N)- gp(N-1)*(1-n1(i,N-1));
        B(N,N-1) = gp(N-1)*n1(i,N);
    
        for j = 2:N-1
            for k = 2:N-1
                if(j==k)
                    B(j,k) = -hp*n1(i-1,j)-hm*n1(i+1,j)-gp(k-1)*(1-n1(i,j-1))-gm(k+1)*(1-n1(i,j+1));
                elseif(j-k == 1)
                    B(j,k) = gp(k)*n1(i,j);
                elseif(j-k == -1)
                    B(j,k) = gm(k)*n1(i,j);
                end
            end
        end
    
        for j = 1:N
            BU(j) = hm*(1-n1(i,j));
            BL(j) = hp*(1-n1(i,j));
        end
    end
    
    F = F+kron(K1,B)+kron(K2,diag(BU,0))+kron(K3,diag(BL,0));
end

G = zeros(M*N);

for i = 1:M
    K1 = zeros(M);
    K2 = zeros(M);
    K3 = zeros(M);
    K1(i,i) = 1;
    hp = h0*exp(-kappa*(M-i));
    hm = h0*exp(-kappa*(i-1));
    for j = 1:M
        if((j-i)==1)
            K2(i,j) = 1;
        elseif((i-j)==1)
            K3(i,j) = 1;
        end
    end
    
    if(i == 1)
        B(1,1)   = -hm*m1(i+1,1) - gm(2)*(1-m1(i,2));
        B(1,2)   = gm(2)*m1(i,1);
        B(N,N)   = -hm*m1(i+1,N)- gp(N-1)*(1-m1(i,N-1));
        B(N,N-1) = gp(N-1)*m1(i,N);
        
        for j = 2:N-1
            for k = 2:N-1
                if(j==k)
                    B(j,k) = -hm*m1(i+1,j)-gp(k-1)*(1-m1(i,j-1))-gm(k+1)*(1-m1(i,j+1));
                elseif(j-k == 1)
                    B(j,k) = gp(k)*m1(i,j);
                elseif(j-k == -1)
                    B(j,k) = gm(k)*m1(i,j);
                end
            end
        end
    
        for j = 1:N
            BU(j) = hm*(1-m1(i,j));
            BL(j) = hp*(1-m1(i,j));
        end
        
    elseif(i == M)
        B(1,1)   = -hp*m1(i-1,1) - gm(2)*(1-m1(i,2));
        B(1,2)   = gm(2)*m1(i,1);
        B(N,N)   = -hp*m1(i-1,N)- gp(N-1)*(1-m1(i,N-1));
        B(N,N-1) = gp(N-1)*m1(i,N);
        
        for j = 2:N-1
            for k = 2:N-1
                if(j==k)
                    B(j,k) = -hp*m1(i-1,j)-gp(k-1)*(1-m1(i,j-1))-gm(k+1)*(1-m1(i,j+1));
                elseif(j-k == 1)
                    B(j,k) = gp(k)*m1(i,j);
                elseif(j-k == -1)
                    B(j,k) = gm(k)*m1(i,j);
                end
            end
        end
    
        for j = 1:N
            BU(j) = hm*(1-m1(i,j));
            BL(j) = hp*(1-m1(i,j));
        end
        
    else
        B(1,1)   = -hp*m1(i-1,1)-hm*m1(i+1,1) - gm(2)*(1-m1(i,2));
        B(1,2)   = gm(2)*m1(i,1);
        B(N,N)   = -hp*m1(i-1,1)-hm*m1(i+1,N)- gp(N-1)*(1-m1(i,N-1));
        B(N,N-1) = gp(N-1)*m1(i,N);
    
        for j = 2:N-1
            for k = 2:N-1
                if(j==k)
                    B(j,k) = -hp*m1(i-1,j)-hm*m1(i+1,j)-gp(k-1)*(1-m1(i,j-1))-gm(k+1)*(1-m1(i,j+1));
                elseif(j-k == 1)
                    B(j,k) = gp(k)*m1(i,j);
                elseif(j-k == -1)
                    B(j,k) = gm(k)*m1(i,j);
                end
            end
        end
    
        for j = 1:N
            BU(j) = hm*(1-m1(i,j));
            BL(j) = hp*(1-m1(i,j));
        end
    end
    
    G = G+kron(K1,B)+kron(K2,diag(BU,0))+kron(K3,diag(BL,0));
end

L = real(eig(F));
J = real(eig(G));


EE(w) = max(L);
MM(w) = min(L);
FF(w) = max(J);
WW(w) = min(J);

end

figure(2)
plot(K, EE,'b',K,FF,'r','LineWidth',5)
set(gca,'fontsize',20)
%axis([0 2 -1.5 1.5])
xlabel('\Lambda')
ylabel('Largest Eigenvalue')

