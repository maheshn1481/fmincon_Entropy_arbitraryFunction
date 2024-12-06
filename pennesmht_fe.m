%% Pennes Solver
function tempqoi = pennesmht_fe(k_t,w_t,rhocp_t,Km,Htime,deltat)
%%% Following are the inputs for the function
% rhocp_t = Tissue Density [kg/m3] * cp_t = Tissue Specific Heat [J/kg/K]
% k_t = Tissue Thermal Conductivity [W/m/K]
% w_t = Blood Perfusion Rate [1/s]
% Qm_t = Metabolic Heat Generation Rate [W/m3]
% H = Applied Magnetic Field [A/m]
% Md = Domain Magnetization Value of MNP [A/m]
% Keff = Magnetic Anistropy Constant of MNP [J/m3]
% deltat = time step in which the H is constant [s]
Qm_t = 0;                   % [W/m3] Tumor Metabolic Heat

%%% Domain Parameters
RT = 0.005;               % [m] Radius of the Tumor
R = 1.1*RT;               % [m] Radius of the computaional domain
vol_t = (4*pi*RT^3)/3;    % [m3] Tumor Volume
dr = RT/50;               % [m] Spatial discretization step size
r = (0:dr:R)';            % [m] Spatial loactions
N = length(r);            % [-] Number of spatial nodes
t_loc = r<=RT;            % [-] Tumor r indices

%%% Time discretization: Implicit method is employed
dt = 1; % time-step for forward runs
t = 0:length(Htime)*deltat;     % Discrete times
TS = length(t);     % Number of time steps

%%% Assigning the PieceWise Stepped H in deltat at respective times at dt resolution
for i = 1: length(Htime)
    startIndex = (i - 1) * (deltat)/dt + 1; % Calculate the starting index for the current iteration
    endIndex = min(startIndex + (deltat)/dt - 1, TS-1); % Calculate the ending index for the current iteration
    HVector(startIndex:endIndex) = Htime(i); % Amplitude vector in time at a resolution of dt assigned from deltat resolution
end

%%% Ambient Conditions
htc = 10;
T_amb = 24;

%%% Blood Properties
rho_b = 1000;           % [kg/m3] Blood Density
cp_b = 3840;            % [J/kg/K] Blood Specific Heat
k_b = 0.50;             % [W/m/K] Blood Thermal Conductivity
T_b = 37;               % [°C] Blood Arterial Temperature
T_bK = T_b + 273.15;    % [K] Blood Arterial Temperature
T_SS = steady_state_temperature(T_b, htc, T_amb, Qm_t, rho_b, cp_b, w_t, k_t, r);
T_initial = T_SS;     % [°C] Initial Condition Temperature
% T_initial = 37;     % [°C] Initial Condition Temperature
alpha_t = k_t/(rhocp_t); % Tissue Thermal diffusivity [m^2/s]
f = 1.63e5;    % [Hz] Magnetic Field Frequency f = 163 kHz

%%% This section has particle magnetic and physical properties for Fe nanoparticles being used in in-vivo experiments 
%%% Magnetic Nanoparticle Physical and Magnetic Parameters MNP Used: Fe3O4
d_mnp = 60e-9;                 % [m] MNP Diameter d = 60 nm
delta_l = 0.1*d_mnp;                 % [m] Liquid Layer Thickness on the Hard Solid MNP delta = 1 nm
dh_mnp = d_mnp+2*delta_l;       % [m] MNP Hydrodynamic Diameter
%%% domain magnetization
% Mass specific magnetization of MNPs used in mice experiments is 150 emu/g
% conversion of emu/g to [A/m]
% 1 emu/g = 1Am2/kg
% Am2/kg to A/m obtained by multiplying with density in kg/m3
rho_mnp = 7874;        % [kg/m3] MNP Density
Md_measured = 150;                  % [emu/g]
Md = Md_measured*rho_mnp;         % [A/m] Domain Magnetization

%%% This section has particle magnetic and physical properties for Fe3O4 nanoparticles just to test
%%% Magnetic Nanoparticle Physical and Magnetic Parameters MNP Used: Fe3O4
% d_mnp = 16e-9;                 % [m] MNP Diameter d = 60 nm
% delta_l = 0.1*d_mnp;                 % [m] Liquid Layer Thickness on the Hard Solid MNP delta = 1 nm
% dh_mnp = d_mnp+2*delta_l;       % [m] MNP Hydrodynamic Diameter
% %%% domain magnetization
% % Mass specific magnetization of MNPs used in mice experiments is 150 emu/g
% % conversion of emu/g to [A/m]
% % 1 emu/g = 1Am2/kg
% % Am2/kg to A/m obtained by multiplying with density in kg/m3
% rho_mnp = 5180;        % [kg/m3] MNP Density
% Md = 446*1000;         % [A/m] Domain Magnetization


%%% MNP Dose and Magnetic Fluid Parameters
mnp_mf_conc = 100; % [ug/ml] micrograms of MNP/milliliter of water
mnp_mf_conc2 = mnp_mf_conc/1000; %[kg/m^3] kg of MNP/m3 of water
mf_vol = 100; %[ul] microliters of nanofluid injected
mf_vol2 = mf_vol*1e-9; % [m3] of nanofluid
mnp_mass_kg = mnp_mf_conc2*mf_vol2;
mnp_vol = mnp_mass_kg/(rho_mnp); % [m3] of MNP
mu_cf = 1e-3;                   % [Pa.s] Viscocity of the Carrier Fluid
mnp_vol_frac = mnp_vol/mf_vol2;  %[-] Volume Fraction of MNP in MF

%%% Constants used in SAR Calculation
mu0 = pi*4e-7;              % [H/m] Permeability of Free Space
kB = 1.38e-23;              % [J/K] Blotzmann Constant
mnp_svol = pi*(d_mnp^3)/6;  % [m3] Solid Volume of Individual Mangetic Nanoparticle
mnp_hvol = pi*(dh_mnp^3)/6; % [m3] Hydrodynamic Volume of Individual Mangetic Nanoparticle
tau0 = 1e-9;                % [s] Attempt Time for Mangetic Moment relaxation
omega = 2*pi*f;             % [rad/s] Angular frequency of applied magnetic field
alpha_CF = 0.55;            % [-] Correction Factor for MNP Heating in MF and Tissue Mediums

%%% SAR Calculations
gamma = Km*mnp_svol/(kB*T_bK);                    % [-] An intermediate parameters used in further calculations
tauB = (3*mu_cf*mnp_hvol)/(kB*T_bK);                % [s] Brownian Relaxation Time of MNP
tauN = tau0*(sqrt(pi)/2)*(1/sqrt(gamma))*exp(gamma) ;   % [s] Neel Relaxation Time of MNP
% if isinf(tauN)
%     error(['tauN is infinity and gamma = ', num2str(gamma)])
% end
One_by_tauE = 1/tauB+1/tauN;                       % [s] Effective Relaxation Time of MNP
tauE = 1/One_by_tauE;

%%% Implicit Finite Difference Scheme using the Cauchy Boundary Condition on the outer boundary
T = zeros(N, TS);
T(:,1) = T_initial; % Initial Condition
% Coefficient Matrix Three-Diagonal Generation
Lower = zeros(N, 1);
Main  = zeros(N, 1);
Upper = zeros(N, 1);
for i = 1:N
    if i ==1 % Symmetry Node
        Main(i) = 1 + 2*alpha_t*dt/(dr^2) + rho_b*cp_b*w_t*dt/(rhocp_t);
        Upper(i) = -2*alpha_t*dt/(dr^2);
    elseif i == N % Fixed Temperature Node
        Lower(i) = -k_t/dr;
        Main(i) = k_t/dr+htc;
    else % Internal Nodes
        Lower(i) = -alpha_t*dt/(dr^2) + alpha_t*dt/(r(i)*dr);
        Main(i) = 1 + 2*alpha_t*dt/(dr^2) + rho_b*cp_b*w_t*dt/(rhocp_t);
        Upper(i) = -alpha_t*dt/(dr^2) - alpha_t*dt/(r(i)*dr);
    end
end
% Solver Time Iterations
for n = 2:TS
    H = HVector(n-1); %
    if H == 0
        SAR = 0;
    else
        zeta = mu0*Md*H*mnp_svol/(kB*T_bK);                 % [-] Langevin parameter
        X_i = mu0*(Md^2)*mnp_vol_frac*mnp_svol/(3*kB*T_bK); % [-] Initial Magnetic Susceptibility
        X_0 = 3*X_i*(coth(zeta)-1/zeta)/zeta;               % [-] Equilibrium Magnetic Susceptibility
        X_L = X_0*omega*tauE/(1+(omega*tauE)^2);            % [-] Loss Component of Magnetic Susceptibility
        P = pi*mu0*f*X_L*H^2;                               % [W/m3] Heat Generation rate of MNPs in MF
        SAR = P/(rho_mnp*mnp_vol_frac);                     % [W/kg] MNP Specific Absorption Rate
        SAR_grams = SAR/1000;                              % [W/g] MNP Specific Absorption Rate
    end
    % %%% Heat Source by MNP in Tumor
    MNP_conc = mnp_mass_kg/vol_t;      % [kg/m3] Concentration of MNP in Tumor, Assuming MNPs are distributed uniformly and confined within the tumor only
    %%% This is used in previous runs
    % MNP_conc = 4; % 4 [mg/cm3] mnp concentration in tumor higher concentration just for testing
    Q_MNP = alpha_CF*MNP_conc*SAR;  % [W/m3] Heat Generation by MNP in Tissue

    %%% Assign Q_MNP to the center region only
    q_mnp = Q_MNP*t_loc;        % MNP heat generation [W/m^3] within the tumor only
    q_mnp_time(:,n-1) = q_mnp;
    Force = zeros(N,1); % Force vector
    Force(1:N-1)  = T(1:N-1,n-1) + (rho_b*cp_b.*w_t*dt*T_b + Qm_t*dt+ q_mnp(1:N-1)*dt)/(rhocp_t);
    Force(N) = htc*T_amb;
    T(:,n) = thomas_algorithm(Lower, Main, Upper, Force); % Spatio-Temporal Temperature
end


% Spatial integration (integrating over radial locations for each time point)
spatial_integral = trapz(r, T, 1);  % Integrate over the radial direction

% Time integration (integrating the result over time)
total_integral = trapz(t, spatial_integral);  % Integrate over time

% T_sum_time = sum(T,2);      % Summation of Temperature over all Times at each r location
% G = trapz(r, T_sum_time);     % Integration of Temperature over the domain in 3D
tempqoi = total_integral;
%%% Plot the results
q_mnp_t = [q_mnp_time(:,1),q_mnp_time];
% [RR,TT] = meshgrid(r,t);
% figure('Position', [40, 40, 800, 600])
% sgtitle('Fe NP Bioheat Solution')
% subplot(2,2,1)
% plot((1:length(HVector))*dt,HVector)
% xlim([0,t(end)])
% xticks([0:t(end)/4:t(end)]);
% xlabel('Time [s]')
% ylabel('H Amplitude, [A/m]');
% title('Magnetic Field');
% 
% subplot(2,2,2)
% surf(RR,TT,q_mnp_t')
% xlabel('Radial Distance [m]')
% ylabel('Time [s]');
% zlabel('Qmnp [W/m^3]');
% title('Spatio-Temporal Heat Source');
% xlim([0,R])
% xticks([0,RT/2,RT,R]);
% ylim([0,t(end)])
% yticks([0:t(end)/4:t(end)]);
% grid on;
% 
% subplot(2,2,3)
% plot_radius = [0,RT,R];
% legend_String = string(plot_radius)+[" m Tumor Center"," m Tumor Edge"," m Outer Boundary"];
% plot_rad_idx = (plot_radius/dr)+1;
% plot(t,T(plot_rad_idx,:), 'LineWidth',2)
% xlabel('Time, t [s]');
% xlim([0,t(end)])
% xticks([0:t(end)/4:t(end)]);
% cylims = ylim;
% ylim ([floor(cylims(1)),ceil(cylims(2))])
% ylabel('Temperature, T [°C]');
% legend(legend_String, Location='best')
% title('Temperature Elevations');
% 
% grid on;
% 
% subplot(2,2,4)
% plot_time = [0,60,120,t(end)];
% legend_String = string(plot_time)+[" s"," s"," s"," s"];
% plot_ind = (plot_time/dt)+1;
% plot(r,T(:,plot_ind),'LineWidth',2)
% xlabel('Radial distance, r [m]');
% xlim([0,R])
% xticks(0:R/5:R);
% cylims = ylim;
% ylim ([floor(cylims(1)),ceil(cylims(2))])
% ylabel('Temperature, T [°C]');
% legend(legend_String, Location='best')
% title('Spatial Temperature profile');
% grid on;
% alltextObjects = findall(gcf,'-property','FontName');
% set(alltextObjects, 'FontName','Times new Roman', 'FontSize',12)
% saveas(gca,['FeNPPBHESolutionForHmax_',num2str(MNP_conc)],'png')

%%% Tri-Diagonal Matrix Alogorithm for Linear System of Equation Solver
    function x = thomas_algorithm(Lower, Main, Upper, Force)
        nn = length(Force);
        c_star = zeros(nn, 1);
        d_star = zeros(nn, 1);
        x = zeros(nn, 1);

        c_star(1) = Upper(1) / Main(1);
        d_star(1) = Force(1) / Main(1);

        for ii = 2:nn
            c_star(ii) = Upper(ii) / (Main(ii) - Lower(ii) * c_star(ii-1));
            d_star(ii) = (Force(ii) - Lower(ii) * d_star(ii-1)) / (Main(ii) - Lower(ii) * c_star(ii-1));
        end

        x(nn) = d_star(nn);
        for ii = nn-1:-1:1
            x(ii) = d_star(ii) - c_star(ii) * x(ii+1);
        end
    end

%%% Steady-State Temperature Profiles without the external heating, to calculate initial condictions for the transient heat transfer model
    function T_SS = steady_state_temperature(T_b, htc, T_amb, Qm_t, rho_b, cp_b, w_t, k_t, r)
        % Inputs:
        % T_b: Blood temperature (°C)
        % htc: Heat transfer coefficient (W/m^2K)
        % T_amb: Ambient temperature (°C)
        % Qm_t: Metabolic heat generation rate (W/m^3)
        % rho_b: Blood density (kg/m^3)
        % cp_b: Blood specific heat capacity (J/kgK)
        % w_t: Blood perfusion rate (1/s)
        % k_t: Tissue thermal conductivity (W/mK)
        % R: Radius of the domain (m)
        % N: Number of radial grid points

        % Discretization parameters
        delta_r = r(2)-r(1);        % Radial step size
        Num_ele = length(r); % Radial positions

        % Initialize system matrix A and right-hand side vector b
        A = zeros(Num_ele, Num_ele);
        b = zeros(Num_ele, 1);

        % Coefficients for blood perfusion and metabolic heat generation
        perfusion_term = w_t * rho_b * cp_b;

        % Fill the matrix A and vector b for internal nodes (i = 2 to N-1)
        for idx = 2:Num_ele-1
            r_i = r(idx);

            % Coefficients for the discretized equation
            coef_T_ip1 = (k_t / delta_r^2) + (k_t / (r_i *  delta_r)); % T_{i+1}
            coef_T_im1 = (k_t / delta_r^2) - (k_t / (r_i *  delta_r)); % T_{i-1}
            coef_T_i = -(2 * k_t / delta_r^2 + perfusion_term);       % T_i

            % Populate the matrix A and vector b
            A(idx, idx+1) = coef_T_ip1;      % T_{i+1}
            A(idx, idx-1) = coef_T_im1;      % T_{i-1}
            A(idx, idx) = coef_T_i;          % T_i
            b(idx) = -perfusion_term * T_b - Qm_t;
        end

        % Boundary condition at r = 0 (Symmetry: T_1 = T_0)
        A(1, 1) = 1;
        A(1, 2) = -1;
        b(1) = 0;

        % Boundary condition at r = R (Convective boundary condition)
        A(Num_ele, Num_ele-1) = k_t / delta_r;
        A(Num_ele, Num_ele) = -(k_t / delta_r + htc);
        b(Num_ele) = -htc * T_amb;

        % Solve the system of linear equations A*T = b
        T_SS = A \ b;

        % Plot the temperature profile
        % figure;
        % plot(r, T_SS, '-o');
        % xlabel('Radial Position (m)');
        % ylabel('Temperature (°C)');
        % title('Temperature Profile in the Tissue');
        % grid on;
    end
end