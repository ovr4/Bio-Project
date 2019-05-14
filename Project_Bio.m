%% This code determines the force require to remove N amount of cell from a surface with R numner of receptors
N = 100 ; %Number of cell/number of repeats
parameters = [0, 1, 2, 3]; % 0 - velcoity,  1 - # receports, 2 - modulus of receptors, 3 - temperature
a_time = [];
a_stress = [];
a_bond = [];
s_bond = [];
a_force =[];
for params = 1:length(parameters)
    multi = [0.25,0.5, 1, 1.5,1.75];
    Average_stress = [];
    Average_breakage_time = [];
    Average_force_bond = [] ;
    Min_force_bond = [];
    for factor=1:length(multi)
        Max_R = 100; % (10^5)Max number of receptors
        E_dot = 1*10^-4; %% 1/s strain_rate changes with velocity 
        E_cell = 2*10^-3;% Pa ( 2 *10^-6 kpa) Modulus constant for spring in cell membrane
        E_receptor = 10*10^-3 ;% Pa (10 *10^-6 kpa) Modulus constant for spring in receptor
        mu_cell = 10^-3;% (Pa*s) Viscosity of cell
        fullx_receptor = 27*10^-9; % m (nm) max extention of receptor
        fullx_cell=5 *10^-9; %max extention of cell
        bond_area = 10^-18;%area of bond feeling the stress - the region connected to substrate

        %% Parameters for MC process
        ko= 4.0; %% 1/s
        kb = 1.38*10^-23;% bottlzman constant
        T= 298 ; %temperature
        A_bond = 10*10^-17; %% area per bond ^-29 82994002

        timestep = 0.001; %seconds 
        velocity = (5*10^-9); % in m/s play with velocity!!!

        %%% Selecting which parameters to change
        if parameters(params) == 0
            velocity = (5*10^-9) * multi(factor);
        end
        if parameters(params) == 1
            Max_R = 100 * multi(factor);
        end

        if parameters(params) == 2
            E_receptor = 10*10^-3 * multi(factor);
        end

        if parameters(params) == 3
            T= 298 * multi(factor);
        end

        %% initialize total stress for cell i
        Cumulative_stress = [];
        t_final = []; % intitiating time required to break bond for each run
        force_bond = []; % initiating average force/bond
        for i =1:N %% repeat it for multiple cells
            Nbonds = Max_R ;% (0.5,1.5) ; %% starting number of receptors bound
            t=0; %% initiaiting time 
            t_bond = 0; %% initialising t_bond
            t_matrix = []; %% initiating t_matrix
            G_matrix = []; %% intitiating G_matrix
            stress_bond = []; %% recording stress on bonds after full extention to calculate force/bond
            Bonds_matrix = []; %% initiating bond matrix
            total_G_receptor = [];
            All_bonds = [ones(Nbonds,1);zeros(Max_R-Nbonds,1)];
            check = 0; % to determine the entry when max time is reached
            %for t=1:time_extention% While loop until Nbonds = 0?
            while Nbonds ~= 0 %% runs until all bonds are broken
                t = t+1;
                %% Changes at every timestep/delta x
                deltax_receptor = t*velocity*timestep; %changes with time, as cell receptor are strethed
                deltax_cell = t*velocity*timestep; %changes with time as cell is strethed 

                %%check length of receptor
                if deltax_receptor >= fullx_receptor
                    deltax_receptor = fullx_receptor;
                end
                %%check length of cell
                if deltax_cell >= fullx_cell
                    deltax_cell = fullx_cell;
                end

                G_receptor = E_receptor*deltax_receptor/fullx_receptor; %stress from receptor - only the deltax changes
                G_cell= E_cell*deltax_cell/fullx_cell; %stress from cell - only the deltax changes
                G_d = mu_cell*E_dot; %dashpot model for cell - strain rate is constant in this scenario
                G_total = Nbonds*(G_receptor)+G_cell+G_d; %% dependes on number of cells bound

                t_matrix = [t_matrix; t*timestep]; %% recoding time t
                G_matrix = [G_matrix; G_total]; %% recording total stress at time t
                Bonds_matrix = [Bonds_matrix; Nbonds]; %% recording total bonds at time t
                %stress_bond = [stress_bond; G_receptor];
                

                if deltax_cell == fullx_cell && deltax_receptor == fullx_receptor %% if the cell is at full extension, you can break the bonds
                    force_to_break = G_receptor*bond_area;
                    total_G_receptor = [total_G_receptor; G_receptor];
                    if check == 0
                        max_G_entry = length(t_matrix);
                        check = 1;
                    end
                    t_bond = t_bond+1;
                    G_cummulative = sum(G_matrix(max_G_entry:end));
                    k_prob = ko*exp((2*fullx_receptor*G_cummulative*A_bond)/(kb*T)); %% ## Ask Blaire chnages when at full extention (and bonds are broken) and the stochaistic bond breaking comences
                    %bondlist = zeros(1,Nbonds);
                    for bond_b=1:length(All_bonds)
                        P = 1-exp(-k_prob*timestep); %% all cells have the same k_prob at the same time (hold it long enough more break)
                        %P=0.9;
                        if P > rand(1) %% bond breaks
                            All_bonds(bond_b) = 0; %% 0 indicates bond was broken, 1 indicates bond is maintained or formed
                            stress_bond = [stress_bond; sum(total_G_receptor)];
                        end
                    end
                    for bond_f=1:length(All_bonds)
                        if 1 < rand(1) %% bond formed
                            All_bonds(bond_f) = 1;

                        end
                    end
                    %%% recording stress on all bonds after fully extended
                    P;
                    %stress_bond = [stress_bond; G_receptor]; 
                end

                    %Nbonds;
                    %sum(bondlist);
                    Nbonds = sum(All_bonds); %% number of bonds left
            end
            Cumulative_stress = [Cumulative_stress; sum(G_matrix)];
            t_final = [t_final; t_matrix(end)];
            %force_bond = max(stress_bond) + mean(stress_bond(2,end));
            %force_bond = [force_bond; max(stress_bond) + mean(stress_bond(2,end))];
            force_bond = stress_bond;
        end



        Average_stress = [Average_stress, mean(Cumulative_stress)];
        Average_breakage_time = [Average_breakage_time, mean(t_final)];
        %Average_force_bond = [Average_force_bond, mean(force_bond)*bond_area];
        Average_force_bond = [Average_force_bond, mean(force_bond)*bond_area];
        Min_force_bond= [Min_force_bond, force_to_break];


%         figure (factor+params*10), clf
%         plot(t_matrix, G_matrix)
%         xlabel('time (seconds)')
%         ylabel('Stress (Pa)')
%         savefig(factor+params*10)
% 
%         figure (factor+1+params*10), clf
%         hist(Cumulative_stress) %% histogram of stess at a given velocity, temperature, and N cells with Max receptors
%         xlabel('Stress (Pa)')
%         savefig(factor+1+params*10)
    end
    a_time = [a_time;Average_breakage_time];
    a_stress = [a_stress; Average_stress] ;
    a_bond = [a_bond; Average_force_bond];
    s_bond = [s_bond; std(Average_force_bond)];
    a_force= [a_force; Min_force_bond];
end
a_time
a_stress
a_bond
a_force
%% Plot properties against velocity

% A_time
figure (1), clf
scatter(multi*velocity*10^9,a_time(1,:), 'filled')
xlabel('velocity (nm/s)')
ylabel('Average breakage Time(s)')

% A-stress
figure (2), clf
scatter(multi*velocity*10^9,a_stress(1,:), 'filled')
xlabel('velocity (nm/s)')
ylabel('Average Stress(Pa)')

%A_bond
parameter_values = [5, 100, 10, 298];
for i = 1:length(parameter_values)
    figure(3*10+i), clf
    scatter(multi*parameter_values(i),a_bond(i,:), 'filled')
    if i == 1
        xlabel('velocity (nm/s)')
    end
    if i == 2
        xlabel('Number of Receptors')
    end
    if i == 3
        xlabel('Modulus (mPa)')
    end
    if i == 4
        xlabel('Temperature (K)')
    end
    ylabel('Average force for bond breaking(N)')

end

%A_Force
figure (4), clf
scatter(multi*velocity*10^9,a_force(1,:), 'filled')
xlabel('velocity (nm/s)')
ylabel('Minimum force for bond(N)')

%% Sensitivity Analysis
parameter_values = [5, 100, 10, 298]; %% /ns, /receptors, /kPa, /K
change= [2,4]; %% 1 is decrease, 3 is increase
a_time_sens=[];
a_stress_sens= [];
a_bond_sens=[];
a_force_sens = [];
for j=1:2
    for i=1:4
        a_time_sens(i,j) = (a_time(i,change(j))-a_time(i,3))/(multi(2)*parameter_values(i));
        a_stress_sens(i,j) = (a_stress(i,change(j))-a_stress(i,3))/(multi(2)*parameter_values(i));
        a_bond_sens(i,j) = (a_bond(i,change(j))-a_bond(i,3))/(multi(2)*parameter_values(i));
        a_force_sens(i,j) = (a_force(i,change(j))-a_force(i,3))/(multi(2)*parameter_values(i));
    end
end

a_time_sens
a_stress_sens
a_bond_sens
a_force_sens

