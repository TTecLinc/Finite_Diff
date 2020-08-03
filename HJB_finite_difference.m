% Peilin Yang
% 2020/8/3

clear all
close all
tic
global w tau ggamma
% setting for saving data and figures
savedata=0; 

% time discount parameters
r=0.05; % annual interest rate
dur_virus   =  1.5; % number of days until cured (cured = vaccine + treatment)
nu = 1/dur_virus ; % daily discount with vaccine + treatment

% propagation parameters
beta    = 0.20 *365 ;  % propagation rate Sdot = beta S .* I / N
ggamma  = 1/18 *365   ;  % avg 20 days to get out of infection state (either recovered or dead)
theta   = 0.5     ;      % effectiveness of transmission    
% parameter for death rate function  Daily Death (t) =  (ggammad +  kappa*I(i)^alphad) .* I(t) ;  
varphi  =  0.01 %   ; % baseline case mortality is 1%  adjusted by duration of desease
kappa   =  0.05 %   ; % benchmark 0.05        % non-linear death rate (units adjusted by duration of desease) increases by 1pp for every 10pp increase in Infected
%%% test tau=1 or no test tau =0
tau =      1    ; % tau =1 test tau=0, no test
%%% max L
l_max = 0.90;  % bench = 0.75;
% Welfare cost parameters
w   =      1    ; % normalized wage
%%% value statistical life
vsl = 20 .* w    ;

%Start 200 grid points
N_S = 200  ;
%N_S=90;
k = 2; % ration of Delta_S = k Delta_I
k=5;
k=5;
N_I = (N_S-1)*k + 1 ;

Delta_S = 1/(N_S-1);
Delta_I = 1/(N_I-1);

S        = linspace(0,1,N_S); 
I        = linspace(0,1,N_I); 
% Add by Peilin
w=2;
h_mi=0;
h_sr=0.2;
b=0.8;
alpha=0.8;
%control=@(S,I)w * (tau *( S+I) + 1-tau );
control=@(S,I,L)(S+ggamma*I)*(log(w*(1/theta)^0.5*L)-0.5+h_sr)+I*(log(b)+h_mi)+vsl*(varphi+kappa*I)*I;
control=@(S,I,L)(S+ggamma*I)*((w+alpha*b)-2*alpha*w*L);
%%% check dt:
max_dt_mat = ones(N_S,N_I)/(r+nu);
for i=2:N_S-1
    for j=2:N_I - (i-1)*k
        max_dt_mat(i,j) = 1/ ( r+nu+ beta* S(i)*I(j).*(1/Delta_S+1/Delta_I) + ggamma* I(j) /Delta_I );
        if S(i)+I(j)>1
            disp('shit')
        end
    end
end

max_dt = min(min(max_dt_mat));

dt = 0.95 .* max_dt

if dt > max_dt
    disp(['dt =  ',num2str(dt)]) ;
    disp('dt is too large');
    disp(['dt must be smaller than ',num2str(max_dt)]);   
end

V = zeros(N_S,N_I);
V_initial = V;
V_next = V ;
L = zeros(N_S,N_I);

for j=1:N_I
     V_initial(1,j) = vsl*I(j) .* ( varphi*ggamma/(r+nu+ggamma) + kappa*ggamma*I(j)/(r+nu+2*ggamma) ) ;
end
for i=2:N_S-1
    for j=2:N_I+1-i
        V_initial(i,j) = control(S(i),I(j),l_max) +  vsl*I(j) .* ( varphi*ggamma/(r+nu+ggamma) + kappa*ggamma*I(j)/(r+nu+2*ggamma) ) ;    
    end
end

%% Test the initial value(Peilin)
% for i=1:N_S
%     for j=1:N_I
%         V_initial(i,j) = w*l_max*( tau*(S(i)+I(j)) +1-tau )/(r+nu) +  vsl*I(j) .* ( varphi*ggamma/(r+nu+ggamma) + kappa*ggamma*I(j)/(r+nu+2*ggamma) ) ;    
%     end
% end

V = V_initial; 
V_next = V_initial; 


diffmax=10; 
diffmax_policy= diffmax;   
iteration=0;
L_optimal = zeros(N_S,N_I);
L_old     = L_optimal ; 
check_foc = zeros(N_S,N_I);
Diff_Der_I_m_S = zeros(N_S,N_I);
Convex_HBJ   = zeros(N_S,N_I); 


%%% Set to NaN to values outside feasible set. 
for i = 1: N_S 
    for j = 1:N_I
        if S(i)+I(j) > 1 %(i+j-2)*Delta>1
            V(i,j)         = NaN;
            V_initial(i,j) = NaN;
            V_next(i,j)   =  NaN;
            Convex_HBJ(i,j) = NaN;
            Diff_Der_I_m_S(i,j) = NaN;
        end
    end
end

N_non_boundary = N_S*N_I- sum(sum(isnan(V_initial))) - N_I - N_S + 1 ;

diff_v = 0; 
s=0;

%while  diffmax>0.0001*r || (diffmax_policy > 0.0001) || iteration < 2/dt 
while iteration<20;
    iteration=iteration+1   ;
    %disp(iteration);
    for i= 2: N_S-1
        %%Step I: Interior of the area;
        for j=2 : N_I-1-(i-1)*k
            %% Step1: get the difference of the grid;
            Diff_Der_I_m_S(i,j)  = 1/Delta_I .* (V(i,j+1) - V(i,j)) - 1/Delta_S .* (V(i,j) - V(i-1,j)) ;
            %disp(Diff_Der_I_m_S(i,j));
            %pause;
            %% Step2.1: right hand side of HJB is convex
            if  Diff_Der_I_m_S(i,j) > 0  
                Convex_HBJ(i,j)= 1 ;
%                 local_A=2*theta^2*Diff_Der_I_m_S(i,j)*beta*S(i)*I(j);
%                 local_B=-2*theta*Diff_Der_I_m_S(i,j)*beta*S(i)*I(j);
%                 local_C=S(i)+I(j)*ggamma;
                Lfoc = ((S(i)+ggamma*I(j))*(w+alpha*b)-2*theta*Diff_Der_I_m_S(i,j)*beta*S(i)*I(j))/...
                     ((S(i)+ggamma*I(j))*(2*w*alpha)-2*theta^2*Diff_Der_I_m_S(i,j)*beta*S(i)*I(j));  %%% use foc
                
                check_foc(i,j) =  control(S(i),I(j),Lfoc) - 2*theta*beta*S(i)*I(j).*(1-theta*Lfoc) .* Diff_Der_I_m_S(i,j) ;
                %L_optimal(i,j) = min(l_max, max(Lfoc, 0));
                L_optimal(i,j)=Lfoc ;
                %L = L_optimal(i,j) ;
                L=Lfoc;
                V_next(i,j)  =   control(S(i),I(j),L) .* dt  +  I(j).*(varphi+kappa*I(j))*ggamma .* vsl .* dt ...
                    + (1-dt*(r+nu)) .* ( 1- ( beta*S(i)*I(j).*(1-theta*L)^2*(dt/Delta_S+dt/Delta_I) + ggamma*I(j).*(dt/Delta_I) ) /(1-dt*(r+nu)) ) .* V(i,j) ...
                    + (1-dt*(r+nu)) .* ( beta*S(i)*I(j).*(1-theta*L)^2 .* (dt/Delta_I)/(1-dt*(r+nu)) ) .* V(i,j+1)   ...
                    + (1-dt*(r+nu)) .* ( beta*S(i)*I(j).*(1-theta*L)^2 .* (dt/Delta_S)/(1-dt*(r+nu)) ) .* V(i-1,j)  ...
                    + (1-dt*(r+nu)) .* ( ggamma*I(j).* (dt/Delta_I)/(1-dt*(r+nu)) ) .* V(i,j-1) ;
                diff_v = diff_v + abs(V_next(i,j)-V(i,j)) ;
                s=s+1;
           
            else 
                %% Step2.2: RHS of HJB is concave
                Convex_HBJ(i,j)= 0 ;
                L = 0;
                V_0  =  control(S(i),I(j),L) .* dt  +  I(j).*(varphi+kappa*I(j))*ggamma .* vsl .* dt ...
                    + (1-dt*(r+nu)) .* ( 1- ( beta*S(i)*I(j).*(1-theta*L)^2*(dt/Delta_S+dt/Delta_I) + ggamma*I(j).*(dt/Delta_I) ) /(1-dt*(r+nu)) ) .* V(i,j) ...
                    + (1-dt*(r+nu)) .* ( beta*S(i)*I(j).*(1-theta*L)^2 .* (dt/Delta_I)/(1-dt*(r+nu)) ) .* V(i,j+1)   ...
                    + (1-dt*(r+nu)) .* ( beta*S(i)*I(j).*(1-theta*L)^2 .* (dt/Delta_S)/(1-dt*(r+nu)) ) .* V(i-1,j)  ...
                    + (1-dt*(r+nu)) .* ( ggamma*I(j).* (dt/Delta_I)/(1-dt*(r+nu)) ) .* V(i,j-1) ;
                L = l_max;
                V_max  =  control(S(i),I(j),L) .* dt  +  I(j).*(varphi+kappa*I(j))*ggamma .* vsl .* dt ...
                    + (1-dt*(r+nu)) .* ( 1- ( beta*S(i)*I(j).*(1-theta*L)^2*(dt/Delta_S+dt/Delta_I) + ggamma*I(j).*(dt/Delta_I) ) /(1-dt*(r+nu)) ) .* V(i,j) ...
                    + (1-dt*(r+nu)) .* ( beta*S(i)*I(j).*(1-theta*L)^2 .* (dt/Delta_I)/(1-dt*(r+nu)) ) .* V(i,j+1)   ...
                    + (1-dt*(r+nu)) .* ( beta*S(i)*I(j).*(1-theta*L)^2 .* (dt/Delta_S)/(1-dt*(r+nu)) ) .* V(i-1,j)  ...
                    + (1-dt*(r+nu)) .* ( ggamma*I(j).* (dt/Delta_I)/(1-dt*(r+nu)) ) .* V(i,j-1) ;
                if V_max < V_0
                    L_optimal(i,j) = l_max;
                    V_next(i,j) = V_max ;
                    diff_v = diff_v + abs(V_next(i,j)-V(i,j)) ;
                    s=s+1;
                else
                    L_optimal(i,j) = 0 ;
                    V_next(i,j) = V_0;
                    diff_v = diff_v + abs(V_next(i,j)-V(i,j)) ;
                    s=s+1;
                end
            end
        end
        %% Step II: Boundary of the area;
        j = N_I - (i-1)*k ;
        Diff_Der_I_m_S(i,j)  = 1/Delta_S .* ( V(i-1,j+k) -V(i,j) ) ;
        %% Step3: optimal L
        if  Diff_Der_I_m_S(i,j) > 0  
            %% Step3.1 right hand side of HJB is convex
            Convex_HBJ(i,j) =1;
%             local_A=2*theta^2*Diff_Der_I_m_S(i,j)*beta*S(i)*I(j);
%             local_B=-2*theta*Diff_Der_I_m_S(i,j)*beta*S(i)*I(j);
%             local_C=S(i)+I(j)*ggamma;
%             Lfoc1 = min((-local_B-(local_B^2-4*local_A*local_C)^0.5/2/local_A),l_max);  %%% use foc
            Lfoc = ((S(i)+ggamma*I(j))*(w+alpha*b)-2*theta*Diff_Der_I_m_S(i,j)*beta*S(i)*I(j))/...
                     ((S(i)+ggamma*I(j))*(2*w*alpha)-2*theta^2*Diff_Der_I_m_S(i,j)*beta*S(i)*I(j));
            %check_foc(i,j) =  control(S(i),I(j),Lfoc) - 2*theta*beta*S(i)*I(j).*(1-theta*Lfoc) .* Diff_Der_I_m_S(i,j) ;
            %L_optimal(i,j) = min(l_max, max(Lfoc, 0));
           L_optimal(i,j)=Lfoc ;
            L=Lfoc;
            V_next(i,j)  =   control(S(i),I(j),L) .* dt  +  I(j).*(varphi+kappa*I(j))*ggamma .* vsl .* dt ...
                + (1-dt*(r+nu)) .* ( 1- ( beta*S(i)*I(j).*(1-theta*L)^2*(dt/Delta_S) + ggamma*I(j)*dt/Delta_I ) /(1-dt*(r+nu)) ) .* V(i,j) ...
                + (1-dt*(r+nu)) .* ( beta*S(i)*I(j).*(1-theta*L)^2 .* (dt/Delta_S)/(1-dt*(r+nu)) ) .* V(i-1,j+k)   ...
                + (1-dt*(r+nu)) .* (  ggamma*I(j) .* (dt/Delta_I)/(1-dt*(r+nu)) )  .* V(i,j-1) ;
            diff_v = diff_v + abs(V_next(i,j)-V(i,j)) ;
            s=s+1;
        else
            %% Step3.2: RHS of HJB is concave
            Convex_HBJ(i,j)= 0 ;
            L = 0;
            V_0  = control(S(i),I(j),L) .* dt  +  I(j).*(varphi+kappa*I(j))*ggamma .* vsl .* dt ...
                + (1-dt*(r+nu)) .* ( 1- ( beta*S(i)*I(j).*(1-theta*L)^2*(dt/Delta_S) + ggamma*I(j)*dt/Delta_I ) /(1-dt*(r+nu)) ) .* V(i,j) ...
                + (1-dt*(r+nu)) .* ( beta*S(i)*I(j).*(1-theta*L)^2 .* (dt/Delta_S)/(1-dt*(r+nu)) ) .* V(i-1,j+k)   ...
                + (1-dt*(r+nu)) .* (  ggamma*I(j) .* (dt/Delta_I)/(1-dt*(r+nu)) )  .* V(i,j-1) ;
            L = l_max;
            V_max  =  control(S(i),I(j),L) .* dt  +  I(j).*(varphi+kappa*I(j))*ggamma .* vsl .* dt ...
                + (1-dt*(r+nu)) .* ( 1- ( beta*S(i)*I(j).*(1-theta*L)^2*(dt/Delta_S) + ggamma*I(j)*dt/Delta_I ) /(1-dt*(r+nu)) ) .* V(i,j) ...
                + (1-dt*(r+nu)) .* ( beta*S(i)*I(j).*(1-theta*L)^2 .* (dt/Delta_S)/(1-dt*(r+nu)) ) .* V(i-1,j+k)   ...
                + (1-dt*(r+nu)) .* (  ggamma*I(j) .* (dt/Delta_I)/(1-dt*(r+nu)) )  .* V(i,j-1) ;
            if V_max < V_0
                L_optimal(i,j) = l_max;
                V_next(i,j) = V_max ;
                diff_v = diff_v + abs(V_next(i,j)-V(i,j)) ;
                s=s+1;
            else
                L_optimal(i,j) = 0 ;
                V_next(i,j) = V_0;
                diff_v = diff_v + abs(V_next(i,j)-V(i,j)) ;
                s=s+1;
            end
        end
    end
    
    diffmax_policy = norm(L_optimal-L_old)/norm(L_old) ;
    diffmax = r*diff_v / N_non_boundary ;
    diff_v = 0;
    % s
    s=0;
    V = V_next ;
    L_old = L_optimal ;
    
end

toc

V_adj = V ; 
V_ini_adj = V_initial;

policy_L = L_optimal;

for i = 1: N_S 
    for j = 1:N_I
        if S(i)+I(j) > 1 %(i+j-2)*Delta>1
            policy_L(i,j)=NaN;
            V_adj(i,j)=NaN;
            V_ini_adj(i,j)=NaN;
        end
    end
end


figure(1)
h2=surf(S,I,policy_L');
set(h2,'Linestyle','none')
th=title('Policy function');set(th,'fontsize',18)
xlabel('Susceptible','FontSize',14); 
ylabel('Infected','FontSize',14); 
zlabel('Policy','FontSize',14); view(0,90)

figure(3)
h2=surf(S,I, V_adj'.*(r)); hold on;
set(h2,'Linestyle','none')
th=title('Value function (flow units)');set(th,'fontsize',18)
xlabel('Susceptible','FontSize',14); 
ylabel('Infected','FontSize',14); 
zlabel('value','FontSize',14); 
view(-55,25)

% figures font size
aaf  = 12 ;   ttf  = 14;  llf  = 14; 
 
figure(4)
subplot(1,2,1)
h2=surf(S,I,policy_L'); hold on;
set(h2,'Linestyle','none')
th=title('Policy function (heat map)');set(th,'fontsize',ttf)
xlabel('Suceptible','FontSize',aaf); 
ylabel('Infected','FontSize',aaf); 
view(0,90)

subplot(1,2,2)
h2=surf(S,I, V_adj'.*(r)); hold on;
set(h2,'Linestyle','none')
th=title('Value function (flow units)');set(th,'fontsize',ttf)
xlabel('Suceptible','FontSize',aaf); 
ylabel('Infected','FontSize',aaf); 
zlabel('value','FontSize',aaf); 
view(-55,25)





