% function [J11,J12,J21,J22] = jacobian_calc(jacobian_params)
%     nbus = jacobian_params.nbus;
% this function does not contain terms multiplied by V in J22 and J12
% meaning Jacobian is not modified, its original one
function [J11,J12,J21,J22] = jacobian_calc_original(jacobian_params)

% this function returns the jacobian when required arguments are passed
    nbus = jacobian_params.nbus;
    G = jacobian_params.G ;
    B = jacobian_params.B;
    Theta = jacobian_params.Theta;
    Y_mag =  jacobian_params.Y_mag;
    PV_bus = jacobian_params.PV_bus;
    Swing_bus = jacobian_params.Swing_bus;
    PQ_bus = jacobian_params.PQ_bus;
    Voltage = jacobian_params.Voltage;
    Delta = jacobian_params.Delta;
% 
% 
% nbus =length(Ybus);
% PV_bus = [4];
% Swing_bus =[1];
% PQ_bus = [2 3];
% Voltage= ones(nbus,1);
% Delta= zeros(nbus,1);
% Ybus =[
% 
%    8.9852-44.8360i -3.8156+19.0781i -5.1696+25.8478i 0.0000+0.0000i;
%   -3.8156+19.0781i 8.9852-44.8360i 0.0000+0.0000i -5.1696+25.8478i;
%   -5.1696+25.8478i 0.0000+0.0000i 8.1933-40.8638i -3.0237+15.1185i;
%    0.0000+0.0000i -5.1696+25.8478i -3.0237+15.1185i 8.1933-40.8638i];
% 
% % getting G and B from Y bus
% G = real(Ybus);
% B = imag(Ybus);
% [Theta Y_mag]=cart2pol(G,B);



    %finding P,Q-calculated
    P_calc = zeros(nbus,1);
    Q_calc = zeros(nbus,1);
    for i=1:nbus
        for j=1:nbus
            P_calc(i) = P_calc(i) + Y_mag(i,j)*Voltage(i)*Voltage(j)*cos(Theta(i,j)+Delta(j)-Delta(i));
            Q_calc(i) = Q_calc(i) - Y_mag(i,j)*Voltage(i)*Voltage(j)*sin(Theta(i,j)+Delta(j)-Delta(i));
        end
    end
    
    
    
    %% calculating J11
    % no need to perform elimination operation later
    nsb = length(Swing_bus); % number of swing bus
    for i=1:nbus
        for j = 1:nbus
            % calculates jacobian only for non-swing buses
            if ~(ismember(i,Swing_bus)| ismember(j,Swing_bus)) 
                if i==j
                    % -nsb is made to correct index of matrix
                    J11(i-nsb,j-nsb) = -Q_calc(i) - (Voltage(i)^2 * B(i,i));
                else
                    J11(i-nsb,j-nsb) = - Y_mag(i,j)*Voltage(i)*Voltage(j)*sin(Theta(i,j)+Delta(j)-Delta(i));
                end
            end
        end
    end
    
    %% calculating J21
    pv_count=0;
    for i=1:nbus
        % if PV bus, dont find jacobian, just increase count required for indexing
        if ismember(i,PV_bus)
            pv_count = pv_count+1;
        else
            for j = 1:nbus
                 % calculated jacobian only for non-swing buses
                if ~(ismember(i,Swing_bus)| ismember(j,Swing_bus));                    
                    if i==j
                        % rows and columns associated with swing bus are eliminated, so nsw is decreased 
                        % from both row and column
                        % apart of nsv, pv_count is subtracted from index as
                        % index gets lowered as well due to exclusion of PV bus
                        %since rows associated to PV of J21 are reduced, index is decreased by pv_count only in row
                        J21(i-nsb-pv_count,j-nsb) = P_calc(i) - (Voltage(i)^2 * G(i,i));
                    else
                        J21(i-nsb-pv_count,j-nsb) = -(Y_mag(i,j)*Voltage(i)*Voltage(j)*cos(Theta(i,j)+Delta(j)-Delta(i)));
                    end
                end
            end
        end
    end
       
    
    
    
    %% calculating J12
    pv_count = 0;
    for j=1:nbus 
        % if PV bus encountered, dont calculate corresponding column elements of J12
        if ismember(j,PV_bus)
                    pv_count=pv_count+1;
        else
            % else calculate element of J12 column wise
            for i = 1:nbus
                if ~(ismember(i,Swing_bus)| ismember(j,Swing_bus)) % if not swing bus
                        if i==j
                            % -nsb, -nsb-pv_count is done in indexing to
                            % arrange indexing

                            %J12(i-nsb,j-nsb-pv_count) = (J21(j-nsb-pv_count,i-nsb) + 2*Voltage(i)^2*G(i,i))/Voltage(i);
                            J12(i-nsb,j-nsb-pv_count) = (P_calc(i) + Voltage(i)^2*G(i,i))/Voltage(i);
                        else
                            %J12(i-nsb,j-nsb-pv_count) = -J21(j-nsb-pv_count,i-nsb)/Voltage(j); mistake
                        
                            J12(i-nsb,j-nsb-pv_count) = Y_mag(i,j)*Voltage(i)*cos(Theta(i,j)+Delta(j)-Delta(i));
                        end
                    end
                end
            end
    end
    
    %% calculating J22
    pv_count_j = 0;
    for j=1:nbus
        if ismember(j,PV_bus)
            % if pv bus is encountered, don't calculate corresponding
            % column of J22
            pv_count_j = pv_count_j + 1;
        else
            pv_count_i=0;
            for i = 1:nbus
                if ~(ismember(i,Swing_bus)| ismember(j,Swing_bus))
                    if ismember(i,PV_bus)
                        % if pv bus is encountered, don't calculate
                        % corresponding column of J22
                        pv_count_i=pv_count_i+1;
                    else
                        if i==j
                            J22(i-nsb-pv_count_i,j-nsb-pv_count_j) = (Q_calc(i) - Voltage(i)^2*B(i,i))/Voltage(i);              
                        else
                          
                            J22(i-nsb-pv_count_i,j-nsb-pv_count_j) = -Y_mag(i,j)*Voltage(i)*sin(Theta(i,j)+Delta(j)-Delta(i));
                        end
                    end
                end
            end
        end
    end
    
    
    J=[J11 J12; J21 J22];
end





