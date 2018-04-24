clear;

%Here we will run a mean field approximation of the spatial moran model for
%a single strain with two possible orientations.  There will be no
%switching between states

%Note that there is a significant amount of commented code at the bottom of
%this page.  Uncommenting this will yield the figure in the supplement.

%N is the number of vertical lattice sites and M is the number of
%horizontal lattice sites.
N = 30;
M = 10;

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

%Strength of the interaction.  Setting to 0 will eliminate preferential growth.
kappa = .012;

%Time step
dt = 0.01;

%Final Time
T  = 8000;

%Creating the hopping matrices.  Since the growth rates or only dependent
%on the row or the column the cell is in based on their orientation, we
%need not have an entire matrix with these quanitities.  Vectors suffice.
for i = 1:M
    hp(i) = h0*exp(-kappa*(M-i));
    hm(i) = h0*exp(-kappa*(i-1));
end

for j = 1:N
    gp(j) = h0*exp(-kappa*(N-j));
    gm(j) = h0*exp(-kappa*(j-1));
end

%Here we use a forward Euler scheme to simulate the full closed system.
for tt = dt:dt:T
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

figure(2)
%subplot(3,3,6)
%colormap('hot')
imagesc(n1)
colorbar
set(gca,'fontsize',20)
axis([1 N 1 M])
xlabel('N')
ylabel('M')
hold on
%zlabel('p_{ij}')
% 
% % %Initialization of the grid
% n1     = rand(M,N);
% n1_new = zeros(M,N);
% %n2 = 1-n1;
% 
% %Setting up the hopping rates.  We note that by taking $\kappa$ = 0, the
% %preferential nature of the hopping will be eliminated
% hp = zeros(M,1);
% hm = zeros(M,1); 
% gp = zeros(N,1);
% gm = zeros(N,1);
% 
% %Basal hopping rate
% h0 = 5;
% 
% %Strength of the interaction.  Setting to 0 will give local sensing model.
% kappa = 0;
% 
% %Time step
% dt = 0.01;
% 
% %Final Time
% T  = 12000;
% 
% %Creating the hopping matrices
% for i = 1:M
%     hp(i) = h0*exp(-kappa*(M-i));
%     hm(i) = h0*exp(-kappa*(i-1));
% end
% 
% for j = 1:N
%     gp(j) = h0*exp(-kappa*(N-j));
%     gm(j) = h0*exp(-kappa*(j-1));
% end
% 
% for tt = dt:dt:T
%     for i = 1:M
%         for j = 1:N
%             if(i ~=1 && i ~= M && j ~= 1 && j~= N)
%                 n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))+hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%             elseif(i == 1)
%                 if(j~=1 && j~=N)
%                    n1_new(i,j) = n1(i,j) + dt*(hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%                 elseif(j==1)
%                    n1_new(i,j) = n1(i,j) + dt*(hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%                 elseif(j==N)
%                    n1_new(i,j) = n1(i,j) + dt*(hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j)));
%                 end
%             elseif(i==M)
%                 if(j~=1 && j~=N)
%                    n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%                 elseif(j==1)
%                    n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%                 elseif(j==N)
%                    n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j)));
%                 end
%             elseif(j==1)
%                 if(i~=1 && i~=M)
%                     n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))+hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%                 end
%             elseif(j==N)
%                 if(i~=1 && i~=M)
%                     n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))+hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j)));
%                 end
%             end
%         end
%     end
%     n1 = n1_new;
% end
% 
% subplot(3,3,4)
% imagesc(n1)
% colorbar
% set(gca,'fontsize',20)
% axis([1 N 1 M])
% xlabel('N')
% ylabel('M')
% 
% 
% %N is the number of vertical lattice sites and M is the number of
% %horizontal lattice sites.
% N = 30;
% M = 10;
% 
% %Initialization of the grid
% n1     = rand(M,N);
% n1_new = zeros(M,N);
% %n2 = 1-n1;
% 
% %Setting up the hopping rates.  We note that by taking $\kappa$ = 0, the
% %preferential nature of the hopping will be eliminated
% hp = zeros(M,1);
% hm = zeros(M,1); 
% gp = zeros(N,1);
% gm = zeros(N,1);
% 
% %Basal hopping rate
% h0 = 5;
% 
% %Strength of the interaction.  Setting to 0 will give local sensing model.
% kappa = .1;
% 
% %Time step
% dt = 0.01;
% 
% %Final Time
% T  = 12000;
% 
% %Creating the hopping matrices
% for i = 1:M
%     hp(i) = h0*exp(-kappa*(M-i));
%     hm(i) = h0*exp(-kappa*(i-1));
% end
% 
% for j = 1:N
%     gp(j) = h0*exp(-kappa*(N-j));
%     gm(j) = h0*exp(-kappa*(j-1));
% end
% 
% for tt = dt:dt:T
%     for i = 1:M
%         for j = 1:N
%             if(i ~=1 && i ~= M && j ~= 1 && j~= N)
%                 n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))+hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%             elseif(i == 1)
%                 if(j~=1 && j~=N)
%                    n1_new(i,j) = n1(i,j) + dt*(hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%                 elseif(j==1)
%                    n1_new(i,j) = n1(i,j) + dt*(hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%                 elseif(j==N)
%                    n1_new(i,j) = n1(i,j) + dt*(hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j)));
%                 end
%             elseif(i==M)
%                 if(j~=1 && j~=N)
%                    n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%                 elseif(j==1)
%                    n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%                 elseif(j==N)
%                    n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j)));
%                 end
%             elseif(j==1)
%                 if(i~=1 && i~=M)
%                     n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))+hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%                 end
%             elseif(j==N)
%                 if(i~=1 && i~=M)
%                     n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))+hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j)));
%                 end
%             end
%         end
%     end
%     n1 = n1_new;
% end
% 
% figure(1)
% subplot(3,3,3)
% %colormap('hot')
% imagesc(n1)
% colorbar
% set(gca,'fontsize',20)
% axis([1 N 1 M])
% xlabel('N')
% ylabel('M')
% hold on
% %zlabel('p_{ij}')
% 
% % %Initialization of the grid
% n1     = rand(M,N);
% n1_new = zeros(M,N);
% %n2 = 1-n1;
% 
% %Setting up the hopping rates.  We note that by taking $\kappa$ = 0, the
% %preferential nature of the hopping will be eliminated
% hp = zeros(M,1);
% hm = zeros(M,1); 
% gp = zeros(N,1);
% gm = zeros(N,1);
% 
% %Basal hopping rate
% h0 = 5;
% 
% %Strength of the interaction.  Setting to 0 will give local sensing model.
% kappa = 0;
% 
% %Time step
% dt = 0.01;
% 
% %Final Time
% T  = 12000;
% 
% %Creating the hopping matrices
% for i = 1:M
%     hp(i) = h0*exp(-kappa*(M-i));
%     hm(i) = h0*exp(-kappa*(i-1));
% end
% 
% for j = 1:N
%     gp(j) = h0*exp(-kappa*(N-j));
%     gm(j) = h0*exp(-kappa*(j-1));
% end
% 
% for tt = dt:dt:T
%     for i = 1:M
%         for j = 1:N
%             if(i ~=1 && i ~= M && j ~= 1 && j~= N)
%                 n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))+hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%             elseif(i == 1)
%                 if(j~=1 && j~=N)
%                    n1_new(i,j) = n1(i,j) + dt*(hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%                 elseif(j==1)
%                    n1_new(i,j) = n1(i,j) + dt*(hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%                 elseif(j==N)
%                    n1_new(i,j) = n1(i,j) + dt*(hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j)));
%                 end
%             elseif(i==M)
%                 if(j~=1 && j~=N)
%                    n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%                 elseif(j==1)
%                    n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%                 elseif(j==N)
%                    n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j)));
%                 end
%             elseif(j==1)
%                 if(i~=1 && i~=M)
%                     n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))+hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%                 end
%             elseif(j==N)
%                 if(i~=1 && i~=M)
%                     n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))+hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j)));
%                 end
%             end
%         end
%     end
%     n1 = n1_new;
% end
% 
% subplot(3,3,1)
% imagesc(n1)
% colorbar
% set(gca,'fontsize',20)
% axis([1 N 1 M])
% xlabel('N')
% ylabel('M')
% 
% %N is the number of vertical lattice sites and M is the number of
% %horizontal lattice sites.
% N = 10;
% M = 30;
% 
% %Initialization of the grid
% n1     = rand(M,N);
% n1_new = zeros(M,N);
% %n2 = 1-n1;
% 
% %Setting up the hopping rates.  We note that by taking $\kappa$ = 0, the
% %preferential nature of the hopping will be eliminated
% hp = zeros(M,1);
% hm = zeros(M,1); 
% gp = zeros(N,1);
% gm = zeros(N,1);
% 
% %Basal hopping rate
% h0 = 5;
% 
% %Strength of the interaction.  Setting to 0 will give local sensing model.
% kappa = .1;
% 
% %Time step
% dt = 0.01;
% 
% %Final Time
% T  = 12000;
% 
% %Creating the hopping matrices
% for i = 1:M
%     hp(i) = h0*exp(-kappa*(M-i));
%     hm(i) = h0*exp(-kappa*(i-1));
% end
% 
% for j = 1:N
%     gp(j) = h0*exp(-kappa*(N-j));
%     gm(j) = h0*exp(-kappa*(j-1));
% end
% 
% for tt = dt:dt:T
%     for i = 1:M
%         for j = 1:N
%             if(i ~=1 && i ~= M && j ~= 1 && j~= N)
%                 n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))+hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%             elseif(i == 1)
%                 if(j~=1 && j~=N)
%                    n1_new(i,j) = n1(i,j) + dt*(hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%                 elseif(j==1)
%                    n1_new(i,j) = n1(i,j) + dt*(hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%                 elseif(j==N)
%                    n1_new(i,j) = n1(i,j) + dt*(hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j)));
%                 end
%             elseif(i==M)
%                 if(j~=1 && j~=N)
%                    n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%                 elseif(j==1)
%                    n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%                 elseif(j==N)
%                    n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j)));
%                 end
%             elseif(j==1)
%                 if(i~=1 && i~=M)
%                     n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))+hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%                 end
%             elseif(j==N)
%                 if(i~=1 && i~=M)
%                     n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))+hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j)));
%                 end
%             end
%         end
%     end
%     n1 = n1_new;
% end
% 
% figure(1)
% subplot(3,3,7)
% %colormap('hot')
% imagesc(n1)
% colorbar
% set(gca,'fontsize',20)
% axis([1 N 1 M])
% xlabel('N')
% ylabel('M')
% hold on
% %zlabel('p_{ij}')
% 
% % %Initialization of the grid
% n1     = rand(M,N);
% n1_new = zeros(M,N);
% %n2 = 1-n1;
% 
% %Setting up the hopping rates.  We note that by taking $\kappa$ = 0, the
% %preferential nature of the hopping will be eliminated
% hp = zeros(M,1);
% hm = zeros(M,1); 
% gp = zeros(N,1);
% gm = zeros(N,1);
% 
% %Basal hopping rate
% h0 = 5;
% 
% %Strength of the interaction.  Setting to 0 will give local sensing model.
% kappa = 0;
% 
% %Time step
% dt = 0.01;
% 
% %Final Time
% T  = 12000;
% 
% %Creating the hopping matrices
% for i = 1:M
%     hp(i) = h0*exp(-kappa*(M-i));
%     hm(i) = h0*exp(-kappa*(i-1));
% end
% 
% for j = 1:N
%     gp(j) = h0*exp(-kappa*(N-j));
%     gm(j) = h0*exp(-kappa*(j-1));
% end
% 
% for tt = dt:dt:T
%     for i = 1:M
%         for j = 1:N
%             if(i ~=1 && i ~= M && j ~= 1 && j~= N)
%                 n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))+hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%             elseif(i == 1)
%                 if(j~=1 && j~=N)
%                    n1_new(i,j) = n1(i,j) + dt*(hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%                 elseif(j==1)
%                    n1_new(i,j) = n1(i,j) + dt*(hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%                 elseif(j==N)
%                    n1_new(i,j) = n1(i,j) + dt*(hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j)));
%                 end
%             elseif(i==M)
%                 if(j~=1 && j~=N)
%                    n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%                 elseif(j==1)
%                    n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%                 elseif(j==N)
%                    n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j)));
%                 end
%             elseif(j==1)
%                 if(i~=1 && i~=M)
%                     n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))+hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%                 end
%             elseif(j==N)
%                 if(i~=1 && i~=M)
%                     n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))+hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j)));
%                 end
%             end
%         end
%     end
%     n1 = n1_new;
% end
% 
% subplot(3,3,9)
% imagesc(n1)
% colorbar
% set(gca,'fontsize',20)
% axis([1 N 1 M])
% xlabel('N')
% ylabel('M')
% 
% %N is the number of vertical lattice sites and M is the number of
% %horizontal lattice sites.
% N = 10;
% M = 30;
% 
% %Initialization of the grid
% n1     = rand(M,N);
% n1_new = zeros(M,N);
% %n2 = 1-n1;
% 
% %Setting up the hopping rates.  We note that by taking $\kappa$ = 0, the
% %preferential nature of the hopping will be eliminated
% hp = zeros(M,1);
% hm = zeros(M,1); 
% gp = zeros(N,1);
% gm = zeros(N,1);
% 
% %Basal hopping rate
% h0 = 5;
% 
% %Strength of the interaction.  Setting to 0 will give local sensing model.
% kappa = .0073;
% 
% %Time step
% dt = 0.01;
% 
% %Final Time
% T  = 12000;
% 
% %Creating the hopping matrices
% for i = 1:M
%     hp(i) = h0*exp(-kappa*(M-i));
%     hm(i) = h0*exp(-kappa*(i-1));
% end
% 
% for j = 1:N
%     gp(j) = h0*exp(-kappa*(N-j));
%     gm(j) = h0*exp(-kappa*(j-1));
% end
% 
% for tt = dt:dt:T
%     for i = 1:M
%         for j = 1:N
%             if(i ~=1 && i ~= M && j ~= 1 && j~= N)
%                 n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))+hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%             elseif(i == 1)
%                 if(j~=1 && j~=N)
%                    n1_new(i,j) = n1(i,j) + dt*(hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%                 elseif(j==1)
%                    n1_new(i,j) = n1(i,j) + dt*(hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%                 elseif(j==N)
%                    n1_new(i,j) = n1(i,j) + dt*(hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j)));
%                 end
%             elseif(i==M)
%                 if(j~=1 && j~=N)
%                    n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%                 elseif(j==1)
%                    n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%                 elseif(j==N)
%                    n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j)));
%                 end
%             elseif(j==1)
%                 if(i~=1 && i~=M)
%                     n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))+hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%                 end
%             elseif(j==N)
%                 if(i~=1 && i~=M)
%                     n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))+hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j)));
%                 end
%             end
%         end
%     end
%     n1 = n1_new;
% end
% 
% figure(1)
% subplot(3,3,8)
% %colormap('hot')
% imagesc(n1)
% colorbar
% set(gca,'fontsize',20)
% axis([1 N 1 M])
% xlabel('N')
% ylabel('M')
% hold on
% %zlabel('p_{ij}')
% 
% %N is the number of vertical lattice sites and M is the number of
% %horizontal lattice sites.
% N = 30;
% M = 10;
% 
% % %Initialization of the grid
% n1     = rand(M,N);
% n1_new = zeros(M,N);
% %n2 = 1-n1;
% 
% %Setting up the hopping rates.  We note that by taking $\kappa$ = 0, the
% %preferential nature of the hopping will be eliminated
% hp = zeros(M,1);
% hm = zeros(M,1); 
% gp = zeros(N,1);
% gm = zeros(N,1);
% 
% %Basal hopping rate
% h0 = 5;
% 
% %Strength of the interaction.  Setting to 0 will give local sensing model.
% %kappa = 0;
% 
% %Time step
% dt = 0.01;
% 
% %Final Time
% T  = 12000;
% 
% %Creating the hopping matrices
% for i = 1:M
%     hp(i) = h0*exp(-kappa*(M-i));
%     hm(i) = h0*exp(-kappa*(i-1));
% end
% 
% for j = 1:N
%     gp(j) = h0*exp(-kappa*(N-j));
%     gm(j) = h0*exp(-kappa*(j-1));
% end
% 
% for tt = dt:dt:T
%     for i = 1:M
%         for j = 1:N
%             if(i ~=1 && i ~= M && j ~= 1 && j~= N)
%                 n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))+hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%             elseif(i == 1)
%                 if(j~=1 && j~=N)
%                    n1_new(i,j) = n1(i,j) + dt*(hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%                 elseif(j==1)
%                    n1_new(i,j) = n1(i,j) + dt*(hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%                 elseif(j==N)
%                    n1_new(i,j) = n1(i,j) + dt*(hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j)));
%                 end
%             elseif(i==M)
%                 if(j~=1 && j~=N)
%                    n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%                 elseif(j==1)
%                    n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%                 elseif(j==N)
%                    n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j)));
%                 end
%             elseif(j==1)
%                 if(i~=1 && i~=M)
%                     n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))+hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gm(j+1)*((1-n1(i,j+1))*n1(i,j)));
%                 end
%             elseif(j==N)
%                 if(i~=1 && i~=M)
%                     n1_new(i,j) = n1(i,j) + dt*(hp(i-1)*n1(i-1,j)*(1-n1(i,j))+hm(i+1)*n1(i+1,j)*(1-n1(i,j))-gp(j-1)*((1-n1(i,j-1))*n1(i,j)));
%                 end
%             end
%         end
%     end
%     n1 = n1_new;
% end
% 
% subplot(3,3,2)
% imagesc(n1)
% colorbar
% set(gca,'fontsize',20)
% axis([1 N 1 M])
% xlabel('N')
% ylabel('M')






