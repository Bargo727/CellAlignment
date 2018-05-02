clear;
%We implement a dynamic stochastic algorithm for cellular dynamics on the lattice.
%We implement a Kinetic Monte-Carlo method based on the Gillespie algorithm.

%System Size
M = 10;
N = 30;

%Parameter values

%Growth Rate
h = 5;

%For kappa > 0, you need to change the mean field formulation.  See
%underneath the full Gillespie loop
kappa = 0;

%Ending Time
FT = 10;
times = 0:.1:FT;
LL = length(times);

%Vector to store the mean fraction of vertical cells
MF = zeros(LL,1);

for ii = 1:LL;
final = times(ii);


%Setting loop
NN = 100000;

%Initiating Gillespie Algorithm

%First, we initialize the lattice to have random initial conditions.  In
%particular, a value of 1 will mean a site is occuupied by a vertical cell
%and a value of 0 indicates occupation by a horizontal cell.
states = [0,1];
L = length(states);

for jj = 1:NN
tt = 0;
    
%Initializing storage vectors

%Vector for fraction of vertical cells
F = [];

%Vector storing times
T = [0];

%Matrix that will store the orientations of the cells at each site.
A = zeros(M,N);

%Random initialization of orientations in lattice.
for j = 1:M
    for i = 1:N
        number = randi([1,L]);
        A(j,i) = states(number);
    end
end

%Initial fraction of vertical cells.  Usually around 0.5.  Will approach
%this value as M,N get large.
F = [F,sum(sum(A))/(M*N)];


while(tt<final)
    
%Vector to store propensities.  Because there are two possible growth
%directions for each cell, we must store them as separate propensities in a
%given vector.  Hence, the vector propensities will be stored in a
%2*M*N-vector. 
p = zeros(M,2*N);
for j = 1:M
    for i = 1:2*N
        if(mod(i,2)==0)
            p(j,i) = (A(j,i/2) == 1)*(j ~= M)*h*exp(-kappa*(M-j)) + (A(j,i/2) == 0)*((i/2)~=N)*h*exp(-kappa*(N-(i/2))); 
        else
            p(j,i) = (A(j,ceil(i/2)) == 1)*(j ~= 1)*h*exp(-kappa*(j-1)) + (A(j,ceil(i/2)) == 0)*(ceil(i/2)~=1)*h*exp(-kappa*(ceil(i/2)-1));
        end
    end
end

P = reshape(p',2*M*N,1);

delta_t = -log(rand(1,1))/sum(P);

%Cumulative sum
selection_intv1 = cumsum(P);

%Normalizing for reaction selection process
selection_intv1 = selection_intv1./selection_intv1(end);

%Selecting index of reaction that will occur
selected_ind = find(selection_intv1 > rand(1,1),1);

%Finding index of corresponding cell
cell_number = ceil(selected_ind/2);

%Finding location of cell that will induce change
c = mod(cell_number,N)*(mod(cell_number,N)~=0) + N*(mod(cell_number,N) ==0);
r = ceil(cell_number/N);

%Parity of selected_ind will determine if growth is toward higher boundary
%or lower boundary.

state = A(r,c) > 0;

if(state > 0)
    if(mod(selected_ind,2)==0)
        A(r+1,c) = A(r,c);
    else
        A(r-1,c) = A(r,c);
    end
else
    if(mod(selected_ind,2)==0)
        A(r,c+1) = A(r,c);
    else
        A(r,c-1) = A(r,c);
    end
    
end

tt = tt + delta_t;

%     if(mod(floor(tt),500)==0)
%         figure(1)
%         imagesc(A)
%         colormap(jet(5))
%         colorbar
%         set(gca,'fontsize',20)
%         drawnow
%     end

F = [F,sum(sum(A))/(M*N)];
T = [T,tt];
end
MF(ii) = MF(ii) + F(end);
end
MF(ii) = MF(ii)/NN;
end

figure(1)
plot(T,F,'r',times,MF,'c','LineWidth',5)
set(gca,'fontsize',20)
xlabel('Time')
ylabel('Fraction of Vertical Cells')
hold on

dt = 0.01;

tt = 0:dt:FT;


%kappa = 0
h1 = h;
g1 = h;

%kappa > 0
%h1 = (h/M)*(1-exp(-kappa*M))/(1-exp(-kappa));
%g1 = (h/N)*(1-exp(-kappa*N))/(1-exp(-kappa));

h0 = h1*(1-(1/(M)));
g0 = g1*(1-(1/(N)));

yy = exp(2*(h0-g0)*tt)./(1 + exp(2*(h0-g0)*tt));

%yy = exp(((1/N)-(1/M))*tt)./(1 + exp(((1/N)-(1/M))*tt));

plot(tt,yy,'k','LineWidth',4)
set(gca,'fontsize',20)

hh = legend('Sample Path','Mean Path','Mean Field');
set(hh,'box','off')

