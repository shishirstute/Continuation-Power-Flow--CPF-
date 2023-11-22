function [P_sch_ori, Q_sch_ori, P_sch, Q_sch] = schedule_power_calc(bus_data, baseMVA,Swing_bus,PV_bus)

    % column 6 and 8 conains real power load and generation
    % this calculates for swing bus also based on bus data but will not be
    % used as P_sch for swing bus is not known before solving
    P_gen = bus_data.data(:,8);
    P_gen = P_gen/baseMVA; % to pu value
    P_load = bus_data.data(:,6);
    P_load = P_load/baseMVA;
    nbus = length(P_load);
    % column 7 and 9 conains reactive power load and generation
    % this calculates for swing bus and PV bus also based on bus data but will not be
    % used as Q_sch for such bus is not known before solving
    Q_gen = bus_data.data(:,9);
    Q_gen = Q_gen/baseMVA; % to pu value
    Q_load = bus_data.data(:,7);
    Q_load = Q_load/baseMVA; % to pu value
    %% finding P,Q-scheduled
    % calculating for P 
    P_sch =  P_gen - P_load;
    % calculating for Q
    Q_sch = Q_gen - Q_load;

    % returning only required P and Q values
    % P and Q for PQ bus, P only for PV bus
   nsb = length(Swing_bus);
   pv_count = 0;
   for i = 1:nbus
       if ~(ismember(i,Swing_bus))
           % for index balancing, nsb is subtracted
           % see jacobian_calc file for further illustration of such
           % concept
           P_sch_ori(i-nsb,1) = P_sch(i);
           if ismember(i,PV_bus)
               pv_count = pv_count +1;
           else
               Q_sch_ori(i-pv_count-nsb,1) = Q_sch(i);
           end
       end
   end
    
end
