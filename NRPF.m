%% iterations starts from here

% maximum iterations is set as 15
% if value converges within prescribed limit, loop will terminate
function [Voltage, Delta, iter] = NRPF(tol_volt,tol_ang,Voltage,Delta,Swing_bus,PQ_bus,PV_bus,nbus...
    ,Y_mag,Theta,bus_data,G,B,baseMVA,bts,lambda,maxiters)

    Non_swing_bus = union(PQ_bus, PV_bus);
    for iter=1:maxiters
    
        if iter > 2
            max_error_volt = max(abs(Voltage_history(:,iter-1)-Voltage_history(:,iter-2)));
            max_error_ang = max(abs(Delta_history(:,iter-1)-Delta_history(:,iter-2)));
            if max_error_volt < tol_volt & max_error_ang < tol_ang
                break;
            end
        end
        % for observing iteration
        iter
        
        % storing v and angle of every iterations
        Voltage_history(:,iter) = Voltage;
        Delta_history(:,iter) = Delta*180/pi; % to degree as well
        
        %% calculating mismatch vector
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
        
        %% calculating jacobian matrix
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
        %[J11,J12,J21,J22] = jacobian_calc(jacobian_params);
        J = [J11 J12;J21 J22];
        
        % solving for delta x
        % delta_x = [delta_angle for non-swing bus; deltaV for PQ bus]
        % passing to solver
        delta_x = crout_solver(J,del_PQ);
        % getting delta_correct(angle) and voltage_correct
        delta_correct = delta_x([1:nbus-length(Swing_bus)]);
        voltage_correct = delta_x(nbus-length(Swing_bus)+1:end);
        
        %% updating
        % creating parameter list
        % update_params.delta_correct = delta_correct;
        % update_params.voltage_correct = voltage_correct;
        % update_params.Swing_bus = Swing_bus;
        % update_params.PV_bus = PV_bus;
        % update_params.Voltage = Voltage;
        % update_params.Delta = Delta;
        % 
        % % calling function
        % [Voltage, Delta] = update_value(update_params);

        % getting del_delta, del_voltage and del_lambda
        Delta(Non_swing_bus)  = Delta(Non_swing_bus) + delta_correct;
        Voltage(PQ_bus)  = Voltage(PQ_bus) + voltage_correct;

    end
end
