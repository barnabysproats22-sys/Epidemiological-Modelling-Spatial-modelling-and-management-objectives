%ODE SIR model code

function [Classes] = ODE_SIR_model(para,ICs,mintime,maxtime)

%Run ODE using ODE45
opts = odeset('RelTol',1e-5);
%Read in ICs
IC = [ICs.S_h' ICs.E_h' ICs.I_h' ICs.R_h' ICs.D_h' ICs.A_h' ICs.P_v' ICs.Q_v' ICs.S_v' ICs.E_v' ICs.I_v'];
IC = reshape(IC,[33,1]);
[t, pop] = ode45(@diff_SIRmodel, [mintime:1: maxtime], IC, opts, para);

%Output of ODE34 is 1 vector, read off into Classes variables, each of
%which is a vector, element i contianing info for farm i.
Classes = struct('S_h',pop(:,[1:3]),'E_h',pop(:,[4:6]),'I_h',pop(:,[7:9]),'R_h',pop(:,[10:12]),'D_h',pop(:,[13:15]),'A_h',pop(:,[16:18]),...
    'P_v',pop(:,[19:21]),'Q_v',pop(:,[22:24]),'S_v',pop(:,[25:27]),'E_v',pop(:,[28:30]),'I_v',pop(:,[31:33]),'t',t);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Diff equations

function dPop = diff_SIRmodel(t,pop,para)

%Put pop into matrix
pop = reshape(pop,[3,11]);

%Now each column is a vector S_h etc, elemts for each farm.

%Assign the population matrix into the classes as vectors
S_h=pop(:,1);
E_h=pop(:,2);
I_h=pop(:,3);
R_h=pop(:,4);
D_h=pop(:,5);
A_h=pop(:,6);

P_v=pop(:,7);
Q_v=pop(:,8);
S_v=pop(:,9);
E_v=pop(:,10);
I_v=pop(:,11);

%Now equation set up is infinitly eaier, write in vector notation

%Define the FOIs. Replaced sum with matirx multiplication (Gives Col Vec)
FOI_h = para.a*para.p_h*(S_h./(S_h+E_h+I_h+R_h)).*para.K*I_v;
FOI_v = para.a*para.p_v*S_v.*(para.K*(I_h./(S_h+E_h+I_h+R_h)));

%Write down the ODE system using vector notation
dS_h = para.mu_h*S_h - FOI_h - para.mu_h*S_h;
dE_h = FOI_h - (para.sigma_h + para.mu_h)*E_h;
dI_h = para.sigma_h*E_h - (para.gamma_h+para.mu_h)*I_h;
dR_h = (1-para.delta)*para.gamma_h*I_h - para.mu_h*R_h;
dD_h = para.delta*para.gamma_h*I_h;
dA_h = para.mu_h*(E_h+I_h+R_h);

N_v = P_v + Q_v + S_v + E_v + I_v;

dP_v = para.B_v.*(ones(para.NumberFarms,1)-para.q*(I_v./N_v)) - para.theta*P_v -...
    P_v.*(P_v + Q_v)./para.C;
dQ_v = (para.B_v).*(para.q*(I_v./N_v)) - para.theta*Q_v -...
    Q_v.*(P_v + Q_v)./para.C;
dS_v = para.theta*P_v - FOI_v - para.mu_v*S_v;
dE_v = FOI_v - (para.sigma_v + para.mu_v)*E_v;
dI_v = para.sigma_v*E_v + para.theta*Q_v - para.mu_v*I_v;

%Reshape the derivatives into a column vector
pop = reshape(pop,[33,1]);
dPop = [dS_h;dE_h;dI_h;dR_h;dD_h;dA_h;dP_v;dQ_v;dS_v;dE_v;dI_v]; %Send back to ODE45 which works with vectors

end

end