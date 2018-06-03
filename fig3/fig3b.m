clear;

%Here we will run a mean field approximation of the spatial moran model for
%a single strain with two possible orientations.  There will be no
%switching between states, i.e. there is no reduction to an Ising model.

%N is the number of vertical lattice sites and M is the number of
%horizontal lattice sites.
N = 30;
M = 15;

%Initialization of the grid
n1     = rand(M,N);
n1_new = zeros(M,N);
%n2 = 1-n1;

%Setting up the hopping rates.  We note that by taking $\kappa$ = 0, the
%preferential nature of the hopping will be eliminated
hp = zeros(M,1);
hm = zeros(M,1); 
gp = zeros(N,1);
gm = zeros(N,1);

%Basal hopping rate
h0 = 5;

%Vector storing interatctions.
dk = .0001;
K = .1:-dk:0;
J = length(K);
B = zeros(1,J);

%Time step
dt = 0.01;

%Final Time
T  = 1500;

loop_count = 100;

for bargo = 1:J
    bargo
    kappa = K(bargo);

%Creating the hopping matrices
for i = 1:M
    hp(i) = h0*exp(-kappa*(M-i));
    hm(i) = h0*exp(-kappa*(i-1));
end

for j = 1:N
    gp(j) = h0*exp(-kappa*(N-j));
    gm(j) = h0*exp(-kappa*(j-1));
end


    %Initialization of the grid
    n1     = rand(M,N);
    n1_new = zeros(M,N);

for k = 1:T
    for i = 1:M
        for j = 1:N
            if(i ~=1 && i ~= M && j ~= 1 && j~= N)
                n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))+hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
            elseif(i == 1)
                if(j~=1 && j~=N)
                   n1_new(i,j) = n1(i,j) + dt*(hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
                elseif(j==1)
                   n1_new(i,j) = n1(i,j) + dt*(hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
                elseif(j==N)
                   n1_new(i,j) = n1(i,j) + dt*(hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j)));
                end
            elseif(i==M)
                if(j~=1 && j~=N)
                   n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
                elseif(j==1)
                   n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
                elseif(j==N)
                   n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j)));
                end
            elseif(j==1)
                if(i~=1 && i~=M)
                    n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))+hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
                end
            elseif(j==N)
                if(i~=1 && i~=M)
                    n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))+hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j)));
                end
            end
        end
    end
    n1 = n1_new;
    
end
B(bargo) = sum(sum(n1))/(M*N);
end

figure(2)
plot(K,B,'r','LineWidth',5)
set(gca,'fontsize',20)
xlabel('\kappa')
ylabel('Fraction of Vertical Cells')
ax = gca;
ax.XTick = [];
ax.YTick = [0 1];

