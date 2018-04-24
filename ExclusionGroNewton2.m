clear;

%Here we find the critical kappa value for distinct interaction kernels.

%Setting seed lattice size
N = 10;
M = 5;

%Vector storing size parameter
S = 1:20;
L = length(S);

%Vectors storing critical kappa values for different interaction kernels.
%Kf is for alpha = 0.1.  Kv is for alpha = 0.5.  Kh is for alpha = 1.  and
%Ke is for exponential kernels.

Kf = zeros(1,L);
Kv = zeros(1,L);
Kh = zeros(1,L);
Ke = zeros(1,L);

for j = 1:L
    s = S(j);
    
    
    
    f1 = 0;
    f2 = 0;
  
     h1 = 0;
     h2 = 0;
     
     s1 = 0;
     s2 = 0;
     
     e1 = 0;
     e2 = 0;
  
    syms k
        for jj = 1:s*M
             f1 = f1 + (1/(1 + k*(s*M-jj)^.1));
             s1 = s1 + (1/(1 + k*(s*M-jj)^.5));
             h1 = h1 + (1/(1 + k*(s*M-jj)));
             e1 = e1 + exp(-k*(s*M-jj));
    
        end
        
        for tt = 1:s*N
             f2 = f2 + (1/(1 + k*(s*N-tt)^.1));
             s2 = s2 + (1/(1+k*(s*N-tt)^.5));
             h2 = h2 + (1/(1 + k*(s*N-tt)));
             e2 = e2 + exp(-k*(s*N-tt));
        end
        
    f = (1/(s*M))*(f1 - 1) - (1/(s*N))*(f2 - 1);
    v = (1/(s*M))*(s1-1) - (1/(s*N))*(s2-1);
    h = (1/(s*M))*(h1 - 1) - (1/(s*N))*(h2 - 1);
    e = (1/(s*M))*(e1 - 1) - (1/(s*N))*(e2 - 1);
     
    assume(k>0);
   
  kkf = max(vpasolve(f));
  kkh = max(vpasolve(h));
  kke = max(vpasolve(e));
  kkv = max(vpasolve(v));
    
  Kf(j) = kkf;
  Kh(j) = kkh;
  Kv(j) = kkv;
  Ke(j) = kke;
end

figure(1)
loglog(S,Kf,'ro','MarkerSize',10,'MarkerFaceColor','r')
loglog(S,Kv,'mo','MarkerSize',10,'MarkerFaceColor','m')
loglog(S,Kh,'bo','MarkerSize',10,'MarkerFaceColor','b')
loglog(S,Ke,'yo','MarkerSize',10,'MarkerFaceColor','y')
set(gca,'fontsize',20)
h = legend('\alpha = .1','\alpha = .5','\alpha = 1','Exponential');
xlabel('log s')
ylabel('\kappa')
ax = gca;
ax.XTick = [0 1];
ax.YTick = [0 -2.5 -5];

