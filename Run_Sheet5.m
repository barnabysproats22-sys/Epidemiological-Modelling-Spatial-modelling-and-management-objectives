%% Parameter Set Up. All In Units Days and Km

%Define color map.
CMap=[0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980;0.9290, 0.6940, 0.1250;0.4940 0.1840 0.5560];

%Inportant parameters
mintime=0;
maxtime = 2e4; %400 Found Heuristically 2e4 for abortions to level off
NumberFarms = 3;
C_Death = 250;  
C_Abort = 50;
C_VC = 500;

%All params that don't change with ICs
para = struct('mu_h',1/(10*365),'p_h',1,'sigma_h',1/4,'gamma_h',1/6,'delta',0.3,...
    'a',1/4,'p_v',0.1,'q',0.05,'m',10,'mu_v',1/16,'theta',1/15,'sigma_v',1/14,'b',0.05,...
    'NumberFarms',NumberFarms);

%Create kernel and add to para
load('locations.mat')
d = Loc2Dist(locations);
K = exp(-para.b*d.^2);
para.('K') = K;

%Add Initial Conditons to para.
para.('N_h0') = [1000,500,500];
para.('N_v0') = para.m*para.N_h0;
para.('C') = para.N_v0';     %Transpose to fit with ODE model

%Use into to define ICs.

ICs = struct('S_h',para.N_h0,'E_h',[0 0 0],'I_h',[1 0 0], 'R_h',[0 0 0],'D_h',[0 0 0],'A_h',[0 0 0],...
    'P_v',para.mu_v*para.N_v0/para.theta,'Q_v',[0 0 0],'S_v',para.N_v0,'E_v',[0 0 0], 'I_v',[0 0 0]);

%Use info to define B_v and add to para
para.('B_v') = ICs.P_v'*(para.theta + para.mu_v/para.theta); %Transpose to fit with ODE model

%% Q2(c)
%Classes.S_h(t,i) gives value for farm i at time t.
[Classes] = ODE_SIR_model(para,ICs,mintime,maxtime);

%Plot Infections by Location and Total
figure(1)
hold on
for i = 1:NumberFarms
    plot(Classes.t,Classes.I_h(:,i),'color',CMap(i,:))
title('Local Infections Against Time','FontSize',18)
ylabel('Infected','FontSize',18)
xlabel('Time (Days)','FontSize',18)
set(gca,'FontSize', 25)
end
plot(Classes.t,sum(Classes.I_h,2),'--k')
hold off
legend('Location 1', 'Location 2', 'Location 3','Total')
axis([0 400 0 Inf])

%Plot Cumulative Deaths and Abortions over time as a subplot
figure(2)
subplot(2,1,1)
hold on
for i = 1:NumberFarms
    plot(Classes.t,Classes.D_h(:,i),'color',CMap(i,:))
title('Cumulative Deaths Over Time','FontSize',18)
ylabel('Deaths','FontSize',18)
xlabel('Time (Days)','FontSize',18)
set(gca,'FontSize', 25)
end
plot(Classes.t,sum(Classes.D_h,2),'--k')
hold off
legend('Location 1', 'Location 2', 'Location 3','Total','location','eastoutside')
axis([0 325 0 Inf]) %Axis chosen for when deaths level off.

subplot(2,1,2)
hold on
for i = 1:NumberFarms
    plot(Classes.t,Classes.A_h(:,i),'color',CMap(i,:))
title('Cumulative Abortions Over Time','FontSize',18)
ylabel('Abortions','FontSize',18)
xlabel('Time (Days)','FontSize',18)
end
plot(Classes.t,sum(Classes.A_h,2),'--k')
hold off
legend('Location 1', 'Location 2', 'Location 3','Total','location','eastoutside')
set(gca,'FontSize', 25)
axis([0 Inf 0 Inf])

%% Q2(d)
%When Infections end for each location
%Need to start looking from 200 else count beginning of infection
for i = 1:NumberFarms
    Duration_Index = find(Classes.E_h([200:end],i)+Classes.I_h([200:end],i)<1,1,'first');
    Duration(i) = Classes.t(Duration_Index+200);
end
T_Duration = max(Duration); %Duration is when the disease dies out in the last farm.

%As Abortion and Death variable are cumlitive use final value in D_h and
%A_h class.
C_Total = C_Death*round(sum(Classes.D_h(end,:))) + C_Abort*round(sum(Classes.A_h(end,:)));

%% Q3(a) Set Up
 Detect_Index = find(sum(Classes.D_h,2)>1,1,'first'); %When disease is detected.
 T_Detect = Classes.t(Detect_Index);

paraS = para;  
paraS.('n_v') = 0.1; %0.4;  %Add n_v to paraS. Change to 0.4 when needed.

%Update C for farm 1 only.
paraS.C(1) = (paraS.N_v0(1)*para.mu_v*paraS.n_v^2)/(para.theta^2*(1-paraS.n_v) + para.mu_v);

%Read in values of NC at 7 days after time detect held in classes to be ICs
ICsS = struct();
fn = fieldnames(Classes);
for k=1:numel(fn)
    ICsS.(fn{k}) = Classes.(fn{k})(Detect_Index+7,:);
end

%For static stratergy.
Classes_Stat = ODE_SIR_model(paraS,ICsS,Detect_Index+7-1,maxtime);

paraPM = para;

%Read in values of NC at 15*7 days after having run VC control.
ICsPM = struct();
fn = fieldnames(Classes_Stat);
for k=1:numel(fn)
    ICsPM.(fn{k}) = Classes_Stat.(fn{k})(15*7,:);
end
 
Classes_PM = ODE_SIR_model(paraPM,ICsPM,Detect_Index + 7 + 15*7-1,maxtime);

%% Q3(a) Figures

%Total Infected Adult Mosquitos Acorss 3 initial stratergies
figure(5)
subplot(1,2,1)
hold on
h=[]; %Set empty so legend always matches number of elements in h.

%Plot the elemnts of each class defined over the time periods required.
h(1)=plot(Classes.t,sum(Classes.I_v(:,:),2),'r'); %NC, just Classes

h(2)=plot(Classes.t([1:Detect_Index+7]),sum(Classes.I_v([1:Detect_Index+7],:),2),'b'); %Static stragergy, plot from Classes and Classes_Stat
plot(Classes_Stat.t([1:maxtime-(Detect_Index+7)]),sum(Classes_Stat.I_v([1:maxtime-(Detect_Index+7)],:),2),'b')

h(3) = plot(Classes.t([1:Detect_Index+7]),sum(Classes.I_v([1:Detect_Index+7],:),2),'m'); %VC->NC Plot from Classes, Classes_Stat and Classes_PM
plot(Classes_Stat.t([1:15*7]),sum(Classes_Stat.I_v([1:15*7],:),2),'m')
plot(Classes_PM.t([1:maxtime-(Detect_Index+7+15*7)])-1,sum(Classes_PM.I_v([1:maxtime-(Detect_Index+7+15*7)],:),2),'m')
n_v = paraS.n_v;
legend(h,{'NC', 'VC\rightarrow VC', 'VC\rightarrow NC'})
title(sprintf('Total Infected Adult Mosquitos n_v = %0.1f',n_v))
ylabel('Infected','FontSize',18)
xlabel('Time (Days)','FontSize',18)
set(gca,'FontSize', 25)
axis([0 500 0 Inf])
hold off

%Total Infected Larvae Acorss 3 initial stratergies
subplot(1,2,2)
hold on
h=[];
h(1)=plot(Classes.t,sum(Classes.Q_v(:,:),2),'r');

h(2)=plot(Classes.t([1:Detect_Index+7]),sum(Classes.Q_v([1:Detect_Index+7],:),2),'b');
plot(Classes_Stat.t([1:maxtime-(Detect_Index+7)]),sum(Classes_Stat.Q_v([1:maxtime-(Detect_Index+7)],:),2),'b')

h(3) = plot(Classes.t([1:Detect_Index+7]),sum(Classes.Q_v([1:Detect_Index+7],:),2),'m');
plot(Classes_Stat.t([1:15*7]),sum(Classes_Stat.Q_v([1:15*7],:),2),'m') 
plot(Classes_PM.t([1:maxtime-(Detect_Index+7+15*7)])-1,sum(Classes_PM.Q_v([1:maxtime-(Detect_Index+7+15*7)],:),2),'m')

legend(h,{'NC', 'VC\rightarrow VC', 'VC\rightarrow NC'})
title(sprintf('Total Infected Larvae n_v = %0.1f',n_v))
ylabel('Infected','FontSize',18)
xlabel('Time (Days)','FontSize',18)
set(gca,'FontSize', 25)
axis([0 500 0 Inf])
hold off

%% Q3(a) Metric Calculation
%Duration NC->VC
for i = 1:NumberFarms
    Duration_IndexS = find(Classes_Stat.E_h([100:end],i)+Classes_Stat.I_h([100:end],i)<1,1,'first');
    DurationS(i) = Classes.t(Duration_IndexS+100+T_Detect+7); %Add Detect +7 to put same time as original classes
end
T_DurationS = max(DurationS);
C_Total_Static = C_Death*(sum(Classes_Stat.D_h(end,:))) + C_Abort*(sum(Classes_Stat.A_h(end,:)))...
    + C_VC*((T_DurationS-T_Detect-7)/7); %Work in Weeks here as C_VC is in weeks.


%Duration NC->VC->NC
for i = 1:NumberFarms
    Duration_Index_PM = find(Classes_PM.E_h([100:end],i)+Classes_PM.I_h([100:end],i)<1,1,'first');
    DurationPM(i) = Classes.t(Duration_Index_PM+100+T_Detect+7+15*7); %Add Detect +7 to put same time as original classes
end
T_DurationPM = max(DurationPM);
C_Total_PM = C_Death*(sum(Classes_PM.D_h(end,:))) + C_Abort*(sum(Classes_PM.A_h(end,:)))...
    +C_VC*15; %Work in Weeks here as C_VC is in weeks.

%% Q3(b)

%n_v varying over all possible values
n_v_vec = [0.05:0.01:1]; %If start at zero get an error. Takes longer to calc closer to zero we get.

%Generate appropriate "f" for integral.
x_0=1;
y_0=4;
f = betapdf(n_v_vec,1+x_0,1+y_0);

%Required if want to alter interval length
T_DurationS_vec =[];
C_Total_Static_vec = [];
T_DurationPM_vec = [];
C_Total_PM_vec = [];
%Generate all Metrics from VC->VC and VC->NCmodel.
for r = 1:length(n_v_vec)
    Detect_Index = find(sum(Classes.D_h,2)>1,1,'first'); %Always same start from Classes.
    T_Detect = Classes.t(Detect_Index);

    paraS = para;  
    paraS.('n_v') = n_v_vec(r);
    paraS.C(1) = (paraS.N_v0(1)*para.mu_v*paraS.n_v^2)/(para.theta^2*(1-paraS.n_v) + para.mu_v);

    %Define the ICs from Classes elvaulated at T_detect+7
    ICsS = struct();
    fn = fieldnames(Classes);
    for k=1:numel(fn)
        ICsS.(fn{k}) = Classes.(fn{k})(Detect_Index+7,:);
    end
    
    ICsS.E_h;
    ICsS.I_h;

    Classes_Stat = ODE_SIR_model(paraS,ICsS,Detect_Index+7-1,maxtime);
   
    for i = 1:NumberFarms
        Duration_IndexS = find(Classes_Stat.E_h([100:end],i)+Classes_Stat.I_h([100:end],i)<1,1,'first');
        DurationS(i) = Classes.t(Duration_IndexS+100+T_Detect+7); %Add Detect +7 to put same time as original classes
    end
    %Store metric values for each n_v.
    T_DurationS_vec(r) = max(DurationS);
    DeathS_vec(r) = sum(Classes_Stat.D_h(end,:));
    AbortS_vec(r) = sum(Classes_Stat.A_h(end,:));
    C_Total_Static_vec(r) = C_Death*(sum(Classes_Stat.D_h(end,:))) + C_Abort*(sum(Classes_Stat.A_h(end,:)))...
        + C_VC*((T_DurationS_vec(r)-T_Detect-7)/7);
    
    %NC->VC->NC Strat
    
    %Define ICs from Classes_Stat after 15*7 days.
    paraPM = para;
    ICsPM = struct();
    fn = fieldnames(Classes_Stat);
    for k=1:numel(fn)
        ICsPM.(fn{k}) = Classes_Stat.(fn{k})(15*7,:);
    end

    Classes_PM = ODE_SIR_model(paraPM,ICsPM,Detect_Index + 7 + 15*7-2,maxtime);

    %Duration NC->VC->NC
    for i = 1:NumberFarms
        Duration_Index_PM = find(Classes_PM.E_h([100:end],i)+Classes_PM.I_h([100:end],i)<1,1,'first');
        DurationPM(i) = Classes.t(Duration_Index_PM+100+T_Detect+7+15*7); %Add Detect +7 to put same time as original classes
    end
    %Store metrics in a vector for each n_v
    T_DurationPM_vec(r) = max(DurationPM);
    DeathPM_vec(r) = sum(Classes_PM.D_h(end,:));
    AbortPM_vec(r) = sum(Classes_PM.A_h(end,:));
    C_Total_PM_vec(r) = C_Death*(sum(Classes_PM.D_h(end,:))) + C_Abort*(sum(Classes_PM.A_h(end,:)))...
        + C_VC*15;

end

%Compute expectations by using trapz to integrate.
Exp_DurationS = round(trapz(n_v_vec,T_DurationS_vec.*f));
Exp_CostS = round(trapz(n_v_vec,C_Total_Static_vec.*f));
Exp_DurationPM = round(trapz(n_v_vec,T_DurationPM_vec.*f));
Exp_CostPM = round(trapz(n_v_vec, C_Total_PM_vec.*f));

%% Q3(c)

%Plot Duration against n_v
figure(7)
clf
subplot(1,2,1)
hold on
plot(n_v_vec, T_DurationS_vec)
plot(n_v_vec, T_DurationPM_vec)
hold off
title('Duration as Varies With n_v','FontSize',18)
xlabel('n_v','FontSize',18)
ylabel('Time (Days)','FontSize',18)
legend('VC\rightarrow VC','VC\rightarrow NC')
set(gca,'FontSize', 25)
axis([-Inf Inf -Inf Inf])

%Plot Cost against n_v
subplot(1,2,2)
hold on
plot(n_v_vec, C_Total_Static_vec/1e3)
plot(n_v_vec, C_Total_PM_vec/1e3)
hold off
title('Cost as Varies With n_v','FontSize',18')
xlabel('n_v','FontSize',18)
ylabel('Cost (1,000 $)','FontSize',18)
legend('VC\rightarrow VC','VC\rightarrow NC')
set(gca,'FontSize', 25)
axis([-Inf Inf -Inf Inf])

%% Plotting the Kernel and Farm Distances

%Plot kernel
figure(8)
subplot(1,2,1)
x = [0:0.01:15];
plot(x,exp(-para.b*x.^2))
title('Plot of Distance Kernel','FontSize',32)
ylabel('K(d)','FontSize',25)
xlabel('d','FontSize',25)
set(gca,'FontSize', 25)

%Plot farm locations
subplot(1,2,2)
hold on
for i = 1:NumberFarms
    plot(locations(i,1),locations(i,2),'.','markersize',40,'color',CMap(i,:));
end
axis([0 20 0 8])
hold off
title('Locations of Farms','FontSize',32)
ylabel('Y (km)','FontSize',25)
xlabel('X (km)','FontSize',25)
set(gca,'FontSize', 25)

%% Plotting C tilde
figure(9)
x = [0:0.001:1];
plot(x,(para.N_v0(1)*para.mu_v*x.^2)./(para.theta^2*(1-x) + para.mu_v))
title('C tilde against n_v','FontSize',32)
ylabel('C tilde','FontSize',25)
xlabel('n_v','FontSize',25)
set(gca,'FontSize', 25)




