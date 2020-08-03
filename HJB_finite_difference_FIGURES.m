% Peilin Yang
% 2020/8/3

clear all
close all

% This code loads an pre-existing output file, runs a time-path simulation and constructs figures 

load('benchmark_vsl_20_gridsize200k5_finitedifferences','-mat');

disp('')
disp([' tau = ',num2str(tau)])
disp([' theta = ',num2str(theta)])
disp([' beta/day = ',num2str(beta/365)])
disp([' gamma/day = ',num2str(ggamma/365)])
disp([' w = ',num2str(w)])
disp([' vsl = ',num2str(vsl)])
disp([' varphi = ',num2str(varphi)])
disp([' kappa = ',num2str(kappa)])
disp([' L max = ',num2str(l_max)])
disp([' r = ',num2str(r)])
disp([' nu = ',num2str(nu)])
disp([' N_I grid points = ',num2str(N_I)])
disp([' N_S grid points = ',num2str(N_S)])
disp([' dt = ',num2str(dt)])
disp([' dt x 365 = ',num2str(dt*365)])
disp('')


figure(1)
h2=surf(S,I,policy_L'); 
set(h2,'Linestyle','none')
th=title('Policy function');set(th,'fontsize',18)
xlabel('Susceptible','FontSize',14); 
ylabel('Infected','FontSize',14); 
zlabel('Policy','FontSize',14); view(0,90)

% value function
figure(3)
h2=surf(S,I, V_adj'*(r)); hold on;
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
h2=surf(S,I, V_adj'*(r)); hold on;
set(h2,'Linestyle','none')
th=title('Value function (flow units)');set(th,'fontsize',ttf)
xlabel('Suceptible','FontSize',aaf); 
ylabel('Infected','FontSize',aaf); 
zlabel('value','FontSize',aaf); 
view(-55,25)

figure(1)
h2=surf(S,I,policy_L');
set(h2,'Linestyle','none')
th=title('Policy function');set(th,'fontsize',18)
xlabel('Susceptible','FontSize',14); 
ylabel('Infected','FontSize',14); 
zlabel('Policy','FontSize',14); view(0,90)

% value function
figure(3)
h2=surf(S,I, V_adj'*(r)); hold on;
set(h2,'Linestyle','none')
th=title('Value function (flow units)');set(th,'fontsize',18)
xlabel('Susceptible','FontSize',14); 
ylabel('Infected','FontSize',14); 
zlabel('value','FontSize',14); 
view(-55,25)

% figures font size
aaf  = 12 ;   ttf  = 14;  llf  = 14; 

%%%%%%%%%%%%%%%%%%

figure(4)
subplot(1,2,1)
h2=surf(S,I,policy_L'); hold on;
set(h2,'Linestyle','none')
th=title('Policy function (heat map)');set(th,'fontsize',ttf)
xlabel('Suceptible','FontSize',aaf); 
ylabel('Infected','FontSize',aaf); 
view(0,90)

subplot(1,2,2)
h2=surf(S,I, V_adj'*(r)); hold on;
set(h2,'Linestyle','none')
th=title('Value function (flow units)');set(th,'fontsize',ttf)
xlabel('Suceptible','FontSize',aaf); 
ylabel('Infected','FontSize',aaf); 
zlabel('value','FontSize',aaf); 
view(-55,25)


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulation with policy function
Nyears = 2 ; % for simulations and figures%

%Dt= dt ;  
Dt = 1/(1*365);
time = round(iteration*dt/Dt) ; 
period=[1:1:time];

% setup vars
Stime=[]; Itime=[]; Rtime=[]; Ntime=[]; Dtime=[];  
%
Stime_Lck=[]; Itime_Lck=[]; Rtime_Lck=[]; Ntime_Lck=[];  Dtime_Lck=[]; 

% initial conditions %%% INITIAL CONDITION
Itime(1)=0.01;
Stime(1)=1-Itime(1)-0.02; % Itime defined above

%%% put it in a grid
[I0,I0_indx] = min(abs(I-Itime(1)));
[S0,S0_indx] = min(abs(S-Stime(1)));

Itime(1) = I(I0_indx);
Stime(1) = S(S0_indx);


Rtime(1)=1-Itime(1)-Stime(1);
Dtime(1)=0;
Ntime(1)=Itime(1)+Rtime(1)+Stime(1);

Itime_Lck(1)=Itime(1); 
Stime_Lck(1)=Stime(1); 
Rtime_Lck(1)=Rtime(1);
Dtime_Lck(1)=0;
Ntime_Lck(1) = Itime_Lck(1)+ Rtime_Lck(1) + Stime_Lck(1);

[I0,I0_indx] = min(abs(I-Itime(1)));
[S0,S0_indx] = min(abs(S-Stime(1)));

Lck(1) = policy_L(S0_indx,I0_indx) ;  

for t=2:time 
    % benchmark uncontrolled process
    Stime(t) =  Stime(t-1)*(1 - Dt*beta*Itime(t-1)*(1-0) );
    Itime(t) =  Itime(t-1)*(1 - Dt*ggamma + Dt*beta*Stime(t-1)*(1-0)) ;

    fatality =  varphi*ggamma + ggamma*kappa*Itime(t-1);
    Dtime(t) =  Dtime(t-1) + Itime(t-1)*fatality*Dt ;
    Rtime(t) =  Rtime(t-1) + Itime(t-1)*Dt*(ggamma-fatality) ;
    Ntime(t) =  Stime(t) + Itime(t) + Rtime(t); 
    
    % Controlled process
    Stime_Lck(t) =  Stime_Lck(t-1)*(1 - Dt*beta*Itime_Lck(t-1)*(1-theta*Lck(t-1))^2) ;
    Itime_Lck(t) =  Itime_Lck(t-1)*(1 - Dt*ggamma + Dt*beta*Stime_Lck(t-1)*(1-theta*Lck(t-1))^2) ;
    
    fatality = varphi*ggamma + ggamma*kappa*Itime_Lck(t-1);
    Dtime_Lck(t) =  Dtime_Lck(t-1) + Itime_Lck(t-1)*fatality*Dt ;
    Rtime_Lck(t) =  Rtime_Lck(t-1) + Itime_Lck(t-1)*Dt*(ggamma-fatality) ;
    Ntime_Lck(t) =  Stime_Lck(t)   + Itime_Lck(t) + Rtime_Lck(t);
    
    % Update for policy
    [It,It_indx] = min(abs(I-Itime_Lck(t) ));
    [St,St_indx] = min(abs(S-Stime_Lck(t) ));
    
    Lck(t) = policy_L(St_indx,It_indx) ;
    
    % Update for no policyf
    [It_nopol,It_nopol_indx] = min(abs(I-Itime(t) ));
    [St_nopol,St_nopol_indx] = min(abs(S-Stime(t) ));     
        
                %%% compute the derivatives as finite differences, two sided when possible
                if St_indx > 1 & St_indx < N_S
                    Der_L_S = (policy_L(St_indx+1,It_indx) - policy_L(St_indx-1,It_indx))/(2*Delta_S);
                elseif St_indx == 1
                    Der_L_S = (policy_L(St_indx+1,It_indx) - policy_L(St_indx,It_indx))/(Delta_S);
                else
                    Der_L_S = (policy_L(St_indx,It_indx) - policy_L(St_indx-1,It_indx))/(Delta_S);
                end
                if It_indx > 1 & It_indx < N_I
                    Der_L_I = (policy_L(St_indx,It_indx+1) - policy_L(St_indx,It_indx-1))/(2*Delta_I);
                elseif It_indx == 1
                    Der_L_I = (policy_L(St_indx,It_indx+1) - policy_L(St_indx,It_indx))/(Delta_I);
                else
                    Der_L_I = (policy_L(St_indx,It_indx) - policy_L(St_indx,It_indx-1))/(Delta_I);
                end
                                
    Lck(t) = policy_L(St_indx,It_indx) + Der_L_S*(Stime_Lck(t)-S(St_indx)) + Der_L_I*(Itime_Lck(t)-I(It_indx)) ;
    Lck(t) = min( l_max, max( Lck(t), 0)) ;

    clear Der_L_S Der_L_I      
end

Frac_Pop_Lck= ( tau*(Stime_Lck +  Itime_Lck) + 1-tau ).*Lck;

figure(7)
plot(period,Stime, 'b' , 'LineWidth', 3); hold on;
plot(period,Stime_Lck, 'r' , 'LineWidth', 3);  hold off;
lh=legend('Susceptible (no  control)','Suceptible (control)');set(lh,'fontsize',18)
xh=xlabel('Time (days)','FontSize',14); 
yh=ylabel('Share of the Population','FontSize',14); 

figure(8)
plot(period,Itime, 'b' , 'LineWidth', 3); hold on;
plot(period,Itime_Lck, 'r' , 'LineWidth', 3);  hold off;
lh=legend('Infected (no  control)','Infected (control)');set(lh,'fontsize',18)
xh=xlabel('Time (days)','FontSize',14); 
yh=ylabel('Share of the Population','FontSize',16); 

figure(9)
plot(period,Rtime, 'b' , 'LineWidth', 3); hold on;
plot(period,Rtime_Lck, 'r' , 'LineWidth', 3);  hold off;
lh=legend('Recovered (no  control)','Recovered (control)');set(lh,'fontsize',18)
xh=xlabel('Time (days)','FontSize',14); 
yh=ylabel('Share of the Population','FontSize',14); 

figure(10)
plot(period,Dtime, 'b' , 'LineWidth', 3); hold on;
plot(period,Dtime_Lck, 'r' , 'LineWidth', 3);  hold off;
lh=legend('Dead (no  control)','Dead (control)');set(lh,'fontsize',18)
xh=xlabel('Time (days)','FontSize',14); 
yh=ylabel('Share of the Population','FontSize',14); 

figure(11)
plot(period,Frac_Pop_Lck, 'k' , 'LineWidth', 3); hold on;
lh=legend('% in lockdown');set(lh,'fontsize',18)
xh=xlabel('Time (days)','FontSize',14); 
yh=ylabel('Share of Population in Lockdown','FontSize',16); 
vline([15 30 45 60 90 120 150 180 365]);

tmax = min(Nyears* round(365/2) , time) ;

% set fontsize figures
aaf  = 12 ;  ttf  = 14; llf  = 14; 

figure(12)

subplot(4,1,1)
plot(period(1:tmax),100*Lck(1:tmax), 'k' , 'LineWidth', 3); hold on;
plot(period(1:tmax),100*Frac_Pop_Lck(1:tmax), 'b' , 'LineWidth', 3);  hold off;
th=title('Lockdown policy');set(th,'fontsize',ttf)
lh=legend('Lockdown rate','% of population '); set(lh,'fontsize',llf)
xh=xlabel('Time','FontSize',aaf); 
yh=ylabel('Population (%)','FontSize',aaf); 
vline([15 30 45 60 90 120 150 180 ]);

subplot(4,1,2)
plot(period(1:tmax),100*Itime(1:tmax), 'b' , 'LineWidth', 3); hold on;
plot(period(1:tmax),100*Itime_Lck(1:tmax), 'r' , 'LineWidth', 3);  hold off;
th=title('Infected (%)');set(th,'fontsize',ttf)
lh=legend('Infected (no  control)','Infected (control)');set(lh,'fontsize',llf)
xh=xlabel('Time','FontSize',aaf); 
yh=ylabel('Population (%)','FontSize',aaf); 

subplot(4,1,3)
plot(period(1:tmax),100*(( Itime_Lck(1:tmax)+Stime_Lck(1:tmax)).*Lck(1:tmax) ), 'b' , 'LineWidth', 3);   
th=title(' GDP losses');set(th,'fontsize',ttf)
lh=legend('Susceptible (control)');set(lh,'fontsize',llf)
xh=xlabel('Time ','FontSize',aaf); 
yh=ylabel('Population (%)','FontSize',aaf); 

subplot(4,1,4)
plot(period(1:tmax),100*Dtime(1:tmax), 'b' , 'LineWidth', 3); hold on;
plot(period(1:tmax),100*Dtime_Lck(1:tmax), 'r' , 'LineWidth', 3);  hold off;
th=title('Dead (%)');set(th,'fontsize',ttf)
lh=legend('Dead (no  control)','Dead (control)');set(lh,'fontsize',llf)
xh=xlabel('Time','FontSize',aaf); 
yh=ylabel('Population (%)','FontSize',aaf); 

rho = exp(-(r+nu)*Dt);

for t=1:time
    rho_vec(1,t)=rho^(t-1);
end
 
Dtime_flow(1)=Dtime_Lck(1,1); Dtime_flow2(1)=Dtime_Lck(1,1);
for t=2:time
    Dtime_flow(1,t)             = Dtime_Lck(t)- Dtime_Lck(t-1);
    Dtime_flow2(1,t)            = Dt*(varphi + kappa*Itime_Lck(t-1))*ggamma*Itime_Lck(t-1);
    Dtime_flow_nopolicy(1,t)    = Dtime(t)-Dtime(t-1);
    Dtime_flow2_nopolicy(1,t)   = Dt*(varphi + kappa*Itime(t-1))*ggamma*Itime(t-1);
end
ratio_Dtime_flow=Dtime_flow./Dtime_flow2;
ratio_Dtime_flow_nopolicy=Dtime_flow_nopolicy./Dtime_flow2_nopolicy;

Itime0 = Itime(1);
V_ini =  V(S0_indx,I0_indx);
V_last    =  V_initial(St_indx,It_indx);

PV_Lck_loss   =      Dt * sum(( tau*(Stime_Lck+Itime_Lck)+1-tau)*w .* Lck .* rho_vec);
PV_Death_loss =      sum(Dtime_flow * vsl .*rho_vec);
PV_loss       =      PV_Lck_loss + PV_Death_loss + rho^time*V_last;

output_flowcost = r * PV_Lck_loss; 
flow_value_fct_opt = r * V_ini
flow_pv_optimal    = r * PV_loss

PV_nopolicy_loss  = 0 ;
PV_nopolicy_Death_loss= sum(Dtime_flow_nopolicy * vsl .*rho_vec);
PV_loss_nopolicy= PV_nopolicy_loss + PV_nopolicy_Death_loss;

V_nopolicy_last      =  V_initial(St_nopol_indx,It_nopol_indx);
flow_pv_nopolicy     = r * (PV_loss_nopolicy + rho^time*V_nopolicy_last)


fake_policy=N_S*ones(time,1);
figure(110)
h2=surf(S,I,policy_L'); 
view(2), shading interp; hold on;
plot3(Stime_Lck(1:time),Itime_Lck(1:time),fake_policy(1:time), '-.c' , 'LineWidth', 4); hold on;
scatter(Stime(1),Itime(1),500,'or','filled'); hold on
plot3(Stime(1:time),Itime(1:time),fake_policy(1:time), '<-r' , 'LineWidth', 4); hold off;
set(h2,'Linestyle','none')
xlabel('Susceptible','FontSize',14); 
ylabel('Infected','FontSize',14);   
zlabel('Policy','FontSize',14); view(0,90)

