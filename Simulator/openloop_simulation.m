clear
clearvars
clc

%--------------------------------------------------------------------------
%% Settings
%--------------------------------------------------------------------------
    timespan = 96; %[hr]
    dt = 20; %[s]
    %----------------------------------------------------------------------
    % Well    	   = [Well1    Well2    Well1     Well2     Well1    Well2]
    % PAR    	   = [    --PI--            --GOR--             --WC--    ]
    % UNCERTAINTY  = [ PI1      PI2      GOR1      GOR2      WC1      WC2 ]
    Theta_nom      = [2.51 0.05 0.2];
    
    %----------------------------------------------------------------------
    % Select the type of input disturbance, w_gc
    % NOTE: 1st column of w_gc correspond to value of w_gc in Sm3/hr
    %       2nd column of w_gc correspond to ratio of step lengths
%     type = 'constant';    w_gc_val = [  40000, 1];
%     type = '1step';       w_gc_val = [  40000, 1;
%                                         20000, 1];
    type = '2step';       w_gc_val = [  40000, 1;
                                        20000, 1;
                                        40000, 1];
    %----------------------------------------------------------------------
    % Specify the percentage distribution of w_gc to each well
    % NOTE: Sum should be <= 100%
    % Well    	= [Well1    Well2]
    w_gc_dist 	= 100; %[%]

%--------------------------------------------------------------------------
%% Confirming the settings
%--------------------------------------------------------------------------
    if sum(w_gc_dist) > 100 % if % distribution of w_gc to each well adds up to more than 100%
        error('Percentage distribution of w_gc to each well cannot add up to more than 100%')
    end

    % Displaying the settings
    fprintf('OPENLOOP SIMULATION SETTINGS \n')
    fprintf('	dt          : %g s                  \n', dt)
    fprintf('	timespan	: %g hrs                \n', timespan)
    fprintf('	w_gc        : %s, [%s] Sm3/hr       \n', type, join(string(w_gc_val(:,1))," "))
    fprintf('	w_ga        : [%s]%%*w_gc           \n', join(string(round(w_gc_dist,1))," "))
    fprintf('	Theta_nom  : [%s]  \n', join(string(Theta_nom)," "))

    % Prompting the user to confirm the settings
%     disp(' ')
%     disp('Please make sure the above settings are okay. Then press any key to start the simulation OR press ctrl + c to stop execution.')
%     pause

%--------------------------------------------------------------------------
%% Creating other necessary data
%--------------------------------------------------------------------------
    t = 0:dt:timespan*60*60; %[s]

    % Well     	= [Well1    Well2]
    m_ga    	= 16422; %[kg]
    m_gt      	= 1669; %[kg]
    m_ot      	= 15229; %[kg]
    m0 = [m_ga, m_gt, m_ot];

    w_gc = make_w_gc(type, w_gc_val, t); %[kg/s]

    w_ga = w_gc_dist/100.*w_gc; %[kg/s]

%--------------------------------------------------------------------------
%% Running openloop simulation
%--------------------------------------------------------------------------
    % Initializing the ODE solver
    update_states = @(m0_k, w_ga_k)update_states(m0_k, w_ga_k, Theta_nom, dt);

    % Preallocation
    m = zeros(length(t),length(m0));    
    m(1,:) = m0;
    Theta_OL = Theta_nom.*ones(length(t),1);

    % For each time step within the simulation timespan
    for k = 1:length(t)-1
        m(k+1,:) = update_states(m(k,:), w_ga(k,:));
    end

    [w_op, w_gop, P_a, P_ainj, rho_l, V_G, rho_m, P_tinj, Y_2, rho_ga, w_ginj, P_wf, w_l, w_o, w_g, P_wh, Y_3, w_gp, w_lp, w_wp] = calculate_output(m, Theta_OL);






%% Plotting the results
%--------------------------------------------------------------------------
    timestamps = t/60/60; %[hrs]
    x_label = 'Time [hr]';
    legend_labels = 'well 1';

    figure
    tiledlayout(2,1,'Padding','compact', 'TileSpacing','compact')
    nexttile;
    plot(timestamps,w_op)
    title('Oil production rate, w_{op}');
    ylim('padded'); xlim('tight')
    %ylims = ylim; x = diff(ylim)*0.1; ylim([ylims(1)-x, ylims(2)+x])
    ylabel('[kg/s]')
    grid on
    xlabel('time [hrs]')

    nexttile
    plot(timestamps,w_gop);
    title('Fluid production rate, w_{gop}');
    xlim([timestamps(1), timestamps(end)]);
%     yline(160,'--r')
    ylabel('[kg/s]');
    xlabel('time [hrs]')
    grid on

 %% saving data to matrix

time_hr = timestamps';
% header = {'Time [hrs]', 'W_ga', 'W_ginj', 'W_g', 'W_l', 'P_a', 'P_ainj', 'P_tinj', 'P_wf', 'P_wh', 'Rho_ga', 'Rho_m', 'Rho_l', 'V_G', 'Y_2', 'Y_3', 'W_op', 'W_wp', 'W_gp', 'W_lp', 'W_gop'};
% matrix = [time_hr, w_ga, w_ginj, w_g, w_l, P_a, P_ainj, P_tinj, P_wf, P_wh, rho_ga, rho_m, rho_l, V_G, Y_2, Y_3, w_op, w_wp, w_gp, w_lp, w_gop]; 


header = {'Time [hrs]', 'W_ga [kg/s]', 'P_wf [bar]', 'P_wh [bar]','W_op [kg/s], W_wp [kg/s], W_gp [kg/s]'};
matrix = [time_hr, w_ga, P_wf, P_wh, w_op, w_wp, w_gp];

% filename = 'well1_const_10pc.csv';
% fid = fopen(filename, 'w');
% fprintf(fid, '%s,', header{1:end});
% fprintf(fid, '%s\n', header{end});
% fclose(fid);
% writematrix(matrix, filename, 'WriteMode', 'append', 'Delimiter', ',');


% plot(tbl, "Var1", "Var2")
