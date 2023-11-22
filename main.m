clc;
clear all;
close all;
% bus to study
bts = 13    ;
%% setting parameters

% tolerance error for voltage and angle for consecutive values
tol_volt = 1e-3;
tol_ang = 1e-3;
tap_include = 1;
maxiters = 25;
%% getting data
[bus_data, branch_data] = data_extract();
baseMVA = 100;

%% Y bus formation
% calling function for ybus calculation
Ybus = y_bus_calculation(bus_data, branch_data, tap_include);
% getting G and B from Y bus
G = real(Ybus);
B = imag(Ybus);
% converting to polar
[Theta Y_mag]=cart2pol(G,B);

%% initializing solution

% finding types of bus
nbus = length(Ybus); % total bus number

% find indexing of PV bus
PV_bus = find(bus_data.data(:,3)==2);

% find indexing of swing bus
Swing_bus =find(bus_data.data(:,3)==3);

% find indexing of PQ bus
PQ_bus = find(bus_data.data(:,3)==0);

% flat start; initialization
% assign 1 to all voltage
Voltage= ones(nbus,1);
% fix voltage of swing bus and PV bus to given value
% column 4 of bus data contains bus voltage, you can also use column 11
Voltage(Swing_bus) = bus_data.data(Swing_bus,4);
Voltage(PV_bus) = bus_data.data(PV_bus,4);

% assign 0 to all bus angles
Delta= zeros(nbus,1);
% fix bus angle of Swing bus to given value
% column 5 of bus data contains angle
Delta(Swing_bus) = bus_data.data(Swing_bus,5) * pi/180;
Non_swing_bus = union(PQ_bus, PV_bus);

% schedule power calculation, For Pv only P, for PQ, both P and Q
% P_sch_ori is P for non-swing bus, Q_schi_ori is Q for PQ bus only
[P_sch_ori, Q_sch_ori, P_sch_all, Q_sch_all] = schedule_power_calc(bus_data, baseMVA,Swing_bus, PV_bus);

% getting K vetor
K = [P_sch_ori; Q_sch_ori];

% continuation parameters 
sigma = 0.1; % for first stage
lambda = 0; % for initialization
Voltage_history1 = [];
Delta_history1 = [];
% getting voltage and angle for starting iteration
[Voltage, Delta,~] = NRPF(tol_volt,tol_ang,Voltage,Delta,Swing_bus,PQ_bus,PV_bus,nbus...
   ,Y_mag,Theta,bus_data,G,B,baseMVA,bts,lambda,maxiters);


%% part I %%%% lambda as continuation parameter

Voltage_collects1 = [];
lambda_collects1 = [];
Voltage_collects1(end+1) = Voltage(bts);
lambda_collects1(end+1) = lambda;
Voltage_history1(:,end+1) = Voltage;
Delta_history1(:,end+1) = Delta;


while (1)
    %%%%%%%%% Predictor %%%%%
    
    % calculating jacobian using last voltage and angle
    %listing parameters for jacobian calculation
    jacobian_params.nbus = nbus;
    jacobian_params.G = G;
    jacobian_params.B = B;
    jacobian_params.Theta = Theta;
    jacobian_params.Y_mag = Y_mag;
    jacobian_params.PV_bus = PV_bus;
    jacobian_params.Swing_bus = Swing_bus;
    jacobian_params.PQ_bus = PQ_bus;
    jacobian_params.Voltage = Voltage;
    jacobian_params.Delta = Delta;

    % calling jacobian function
    [J11,J12,J21,J22] = jacobian_calc_original(jacobian_params);
    J = [J11 J12;J21 J22];
    % finding ek
    ek = [zeros(1,length(J)) 1];
    % finding augmented jacobian
    J_aug = [J K;ek]; 

    % getting del_delta, del_voltage and del_lambda using crout based solver
    del_x = crout_solver(J_aug,ek');

    % getting del_delta, del_voltage and del_lambda
    del_delta = del_x([1:nbus-length(Swing_bus)]);
    del_voltage = del_x(nbus-length(Swing_bus)+1:end-1);
    del_lambda = del_x(end);

    % updating
    Delta(Non_swing_bus)  = Delta(Non_swing_bus) + sigma* del_delta;
    Voltage(PQ_bus)  = Voltage(PQ_bus) + sigma* del_voltage;
    lambda = lambda + sigma*del_lambda;

    %%%%%% Corrector %%%%%%%%
    
    % updating studies load value for NR
    [Voltage, Delta, iter] = NRPF(tol_volt,tol_ang,Voltage,Delta,Swing_bus,PQ_bus,PV_bus,nbus...
       ,Y_mag,Theta,bus_data,G,B,baseMVA,bts,lambda,maxiters);
    
    if (iter == maxiters)
        % assign previous value of V and lambda as the solution diverges
        Voltage = Voltage_history1(:,end);
        lambda = lambda_collects1(end);
        Voltage_collects1(end+1) = Voltage(bts);
        lambda_collects1(end+1) = lambda;
        Delta = Delta_history1(:,end);
        break;
    else
        Voltage_history1(:,end+1) = Voltage;
        Delta_history1(:,end+1) = Delta;
        Voltage_collects1(end+1) = Voltage(bts);
        lambda_collects1(end+1) = lambda;
    end
end

%% Part II %%%% Voltage as continuation parameter

sigma = 0.005; % for 12 and 13 use 0.0005
factor = 0.75; % factor where the continuation parameter switches to lambda for 12 and 13 use 0.75
lambda = lambda_collects1(end);
Voltage = Voltage_history1(:,end);
Delta = Delta_history1(:,end);
Voltage_history2 = [];
Delta_history2 = [];
lambda_collects2 = [lambda_collects1(end)];
Voltage_collects2 = [Voltage_collects1(end)];


%%% Predictor step %%%
lambda_2start = lambda;
iteration = 1;
while (lambda>lambda_2start*factor) && iteration < 100
    iteration = iteration+1;
    %finding jacobian
    %listing parameters for jacobian calculation
    jacobian_params.nbus = nbus;
    jacobian_params.G = G;
    jacobian_params.B = B;
    jacobian_params.Theta = Theta;
    jacobian_params.Y_mag = Y_mag;
    jacobian_params.PV_bus = PV_bus;
    jacobian_params.Swing_bus = Swing_bus;
    jacobian_params.PQ_bus = PQ_bus;
    jacobian_params.Voltage = Voltage;
    jacobian_params.Delta = Delta;

    % calling jacobian function
    [J11,J12,J21,J22] = jacobian_calc_original(jacobian_params);
    J = [J11 J12;J21 J22];

    % finding ek
    ek = [zeros(1,length(J)) 1];

    % finding ekv i.e. last row vector of augmented jacobian for part II
    ekv = [zeros(1,length(J)) 0];

    % finding index of voltage continuation variable in augmented jacobian
    modified_Qindex = find(PQ_bus(:)==bts);
    v_index = length(Non_swing_bus) + modified_Qindex;
    ekv(v_index) = -1;
    
    % augmented jacobian
    J_aug = [J K;ekv]; 
    
    % getting del_delta, del_voltage and del_lambda using crout based solver
    del_x = crout_solver(J_aug,ek');
    % getting del_delta, del_voltage and del_lambda
    del_delta = del_x([1:nbus-length(Swing_bus)]);
    del_voltage = del_x(nbus-length(Swing_bus)+1:end-1);
    del_lambda = del_x(end);
    
    % updating
    Delta(Non_swing_bus)  = Delta(Non_swing_bus) + sigma* del_delta;
    Voltage(PQ_bus)  = Voltage(PQ_bus) + sigma* del_voltage;
    lambda = lambda + sigma*del_lambda;
   
    %%%%     Corrector Step %%%%%%%%%%%%%%%%%%%%%%

    jacobian_params.nbus = nbus;
    jacobian_params.G = G;
    jacobian_params.B = B;
    jacobian_params.Theta = Theta;
    jacobian_params.Y_mag = Y_mag;
    jacobian_params.PV_bus = PV_bus;
    jacobian_params.Swing_bus = Swing_bus;
    jacobian_params.PQ_bus = PQ_bus;
    jacobian_params.Voltage = Voltage;
    jacobian_params.Delta = Delta;

    J_aug = [J -K; ekv];

    %  mismatch power calculations
    % listing parameters for mismatch calculation
    mismatch_calc_params.Swing_bus = Swing_bus;
    mismatch_calc_params.PQ_bus = PQ_bus;
    mismatch_calc_params.PV_bus = PV_bus;
    mismatch_calc_params.nbus = nbus;
    mismatch_calc_params.Y_mag = Y_mag;
    mismatch_calc_params.Theta = Theta;
    mismatch_calc_params.Delta = Delta;
    mismatch_calc_params.Voltage = Voltage;
    mismatch_calc_params.bus_data = bus_data;
    mismatch_calc_params.baseMVA = baseMVA;
    mismatch_calc_params.bts = bts;
    mismatch_calc_params.lambda = lambda;

    % calling function
    [del_P del_Q P_calc Q_calc] = mismatch_calc(mismatch_calc_params);
    del_PQ = [del_P; del_Q];

    % getting rhs for corrector
    del_PQV = [del_PQ; 0];

    % passing to solver
    % getting del_delta, del_voltage and del_lambda using crout based solver
    del_x = crout_solver(J_aug,del_PQV);
    % getting del_delta, del_voltage and del_lambda
    del_delta = del_x([1:nbus-length(Swing_bus)]);
    del_voltage = del_x(nbus-length(Swing_bus)+1:end-1);
    del_lambda = del_x(end);
    
    Delta(Non_swing_bus)  = Delta(Non_swing_bus) + del_delta;
    Voltage(PQ_bus)  = Voltage(PQ_bus) + del_voltage;
    lambda = lambda + del_lambda;

    % listings values
    %Voltage = Voltage_pred; % added later for debug
    Voltage_history2(:,end+1) = Voltage;
    Delta_history2(:,end+1) = Delta;
    Voltage_collects2(end+1) = Voltage(bts);
    lambda_collects2(end+1) = lambda;
end
%% Part III Lambda as continuation parameter %%%%

sigma = 0.1;
lambda_collects3 = [lambda_collects2(end)];
Voltage_collects3 = [Voltage_collects2(end)];
Voltage_history3 = [];
Delta_history3 = [];
iteration = 1;
while iteration < 100
    %%%%%%%%% Predictor %%%%%
    iteration = iteration +1;
    % calculating jacobian using last voltage and angle
    %listing parameters for jacobian calculation
    jacobian_params.nbus = nbus;
    jacobian_params.G = G;
    jacobian_params.B = B;
    jacobian_params.Theta = Theta;
    jacobian_params.Y_mag = Y_mag;
    jacobian_params.PV_bus = PV_bus;
    jacobian_params.Swing_bus = Swing_bus;
    jacobian_params.PQ_bus = PQ_bus;
    jacobian_params.Voltage = Voltage;
    jacobian_params.Delta = Delta;
    % calling jacobian function
    [J11,J12,J21,J22] = jacobian_calc_original(jacobian_params);
    J = [J11 J12;J21 J22];
    
    % finding ek
    ek = [zeros(1,length(J)) -1];
    % finding augmented jacobian
    J_aug = [J -K;ek]; 

    % getting del_delta, del_voltage and del_lambda using crout based solver
    del_x = crout_solver(J_aug,abs(ek)');
    % getting del_delta, del_voltage and del_lambda
    del_delta = del_x([1:nbus-length(Swing_bus)]);
    del_voltage = del_x(nbus-length(Swing_bus)+1:end-1);
    del_lambda = del_x(end);
    
    % updating
    Delta(Non_swing_bus)  = Delta(Non_swing_bus) + sigma* del_delta;
    Voltage(PQ_bus)  = Voltage(PQ_bus) + sigma* del_voltage;
    lambda = lambda + sigma*del_lambda;

 
    %% corrector %%%
    [Voltage, Delta, iter] = NRPF(tol_volt,tol_ang,Voltage,Delta,Swing_bus,PQ_bus,PV_bus,nbus...
       ,Y_mag,Theta,bus_data,G,B,baseMVA,bts,lambda,maxiters);
    
    if lambda < 0
        Voltage_collects3(end+1) = Voltage(bts);
        lambda_collects3(end+1) = 0;
        break;
    else
        Voltage_history3(:,end+1) = Voltage;
        Delta_history3(:,end+1) = Delta;
        Voltage_collects3(end+1) = Voltage(bts);
        lambda_collects3(end+1) = lambda;
    end
end


%% plotting
figure;
plot(lambda_collects1, Voltage_collects1, '-r',lambda_collects2,Voltage_collects2,'-g',lambda_collects3,Voltage_collects3,'-k', 'LineWidth', 2.5);
grid on;
xlabel('lambda','Fontsize', 12);
ylabel('Voltage in pu', 'Fontsize',12);
title(['Continuation power flow for bus ' num2str(bts)]);

% folder to save result
folderPath = "C:\Users\shish\OneDrive - Washington State University (email.wsu.edu)\WSU_class\fall 2023\EE521\Homework\Homework2\results";
fileName = ['bus', num2str(bts),'.png'];
fullfilepath = fullfile(folderPath,fileName);
saveas(gcf,fullfilepath);













