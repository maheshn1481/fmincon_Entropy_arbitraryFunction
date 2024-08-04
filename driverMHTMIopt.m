clear all
clc

myoptions.Algorithm = 'constDirect'
%myoptions = optimoptions(@fmincon,'Display','iter-detailed','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'MaxFunctionEvaluations',1e7,'ConstraintTolerance',2.e-9, 'OptimalityTolerance',2.5e-4,'Algorithm','interior-point','StepTolerance',1.000000e-12,'MaxIterations',1000,'PlotFcn',{'optimplotfvalconstr', 'optimplotconstrviolation', 'optimplotfirstorderopt' },'HonorBounds',true, 'Diagnostic','on','FunValCheck','on' )

%% % monitor memory: while [ -e /proc/3291925 ] ; do  top -b -n 1 -p 3291925 >>process.txt ;sleep 60; done

%function driverHPMIopt(NGauss,NumberUncertain,modelSNR,myoptions,ObjectiveType,GaussLegendre )
NGauss = 3,NumberUncertain=1,modelSNR=10, ObjectiveType = 'TotalSignal',GaussLegendre=false


NGauss,NumberUncertain,modelSNR,myoptions.Algorithm,ObjectiveType,GaussLegendre
close all
%% Tissue Parameters
kmean = [ .5 ]; % s
kstdd = [ .2 ]; % s
klb   = [ 5  ]; % s
kub   = [ 45 ]; % s
wmean = [ 6  ]; % s
wstdd = [ 3  ]; % s
wlb   = [ 2  ]; % s
wub   = [ 10 ]; % s
crhomean = [ 30000 ];       % s
crhostdd = [ 10000 ];       % s
crholb   = [ 20000 ];       % s
crhoub   = [ 40000 ];       % s
chimean = [ 0.05 ];      % s
chistdd = [ .01  ];      % s
chilb   = [ 0.01 ];      % s
chiub   = [ 0.20 ];      % s
ubmean  = [ 4    ];      % s
ubsttd  = [ 1.3  ];      % s
ublb    = [ 0    ];       % s
ubub    = [ 7    ];       % s

%% Variable Setup
Ntime = 30;
currentTR = 3;
timelist = (0:(Ntime-1))*currentTR ;
Power = 20*ones(Ntime,1);
temperature = zeros(Ntime,1);

for iii = 1:Ntime
    temperature(iii) = pennesmht(kmean,wmean,crhomean,.001,Power(iii));
end

%% Plot initial guess
plotinit = true;
if plotinit
    % plot initial guess
    figure(1)
    plot(timelist ,temperature)
    ylabel('temperature ')
    xlabel('sec')
end


%% optimize MI
optf = true;
if optf

    tic;
    % setup optimization variables
    Nspecies = 1
    powerList = optimvar('powerList',Nspecies,Ntime,'LowerBound',0, 'UpperBound',35);

    signu = sqrt(2* Ntime) * 10;
    [x2,xn2,xm2,w2,wn2]=GaussHermiteNDGauss(NGauss,0,signu);
    lqp2=length(xn2{1}(:));

    % switch between uniform and Gaussian RV for prior
    if(GaussLegendre)
        QuadratureRule = 'Legendre';
    else
        QuadratureRule = 'Hermite';
    end
    switch (NumberUncertain)
        case(1)
            if(GaussLegendre)
                [x,xn,xm,w,wn]=GaussLegendreNDGauss(NGauss,chilb,chiub);
            else
                [x,xn,xm,w,wn]=GaussHermiteNDGauss(NGauss,chimean,chistdd);
            end
            kqp     = kmean;
            wqp     = wmean;
            crhoqp  = crhomean;
            chiqp   = xn{1}(:);
        case(4)
            if(GaussLegendre)
                [x,xn,xm,w,wn]=GaussLegendreNDGauss(NGauss,[klb; wlb; crholb; chilb],[kub; wub; crhoub; chiub]);
            else
                [x,xn,xm,w,wn]=GaussHermiteNDGauss(NGauss,[kmean; wmean; crhomean; chimean],[kstdd; wstdd; crhostdd; chistdd]);
            end
            kqp    = xn{1}(:);
            wqp    = xn{2}(:);
            crhoqp = xn{3}(:);
            chiqp  = xn{4}(:);
    end
    %alphaqp = xn{6}(:);
    %betaqp  = xn{7}(:);

    lqp=length(xn{1}(:));
    %statevariable    = optimvar('state',Ntime,Nspecies,lqp,'LowerBound',0,'UpperBound',.1);
    statevariable  = optimvar('state',Ntime,Nspecies,lqp,'LowerBound',0);
    auxvariable      = optimexpr(   [Ntime,Nspecies,lqp]);
    stateconstraint  = optimconstr(    [Ntime,Nspecies,lqp]);

    disp('build state variable')

    % [expATRoneone(end),0; expATRtwoone(end), expATRtwotwo(end)]
    % IC
    stateconstraint(1,:,:)  = statevariable(1,:,:) ==0;
    auxvariable(1,:,:) =0;
    for iii = 1:Ntime-1
        % setup state as linear constraint
        auxvariable(iii+1,1,:) =  pennesmht(kqp,wqp,crhoqp,chiqp,powerList(iii+1));
        stateconstraint(iii+1,1,:)  = statevariable(iii+1,1,:) ==  pennesmht(kqp,wqp,crhoqp,chiqp,powerList(iii+1));
    end


    disp('build objective function')
    expandvar  = ones(1,lqp);

    switch (ObjectiveType)
        case('TotalSignal')
            % TODO - repmat does not work well with AD
            % TODO - replace repmat with matrix
            % sumstatevariable = squeeze(sum(repmat(sin(FaList)',1,1,lqp).*statevariable,1));
            sumstatevariable = optimexpr([Nspecies,lqp]);
            for jjj = 1:lqp
                sumstatevariable(:,jjj) =  sum(statevariable(:,:,jjj) ,1)';

            end
            diffsumm =sumstatevariable(1,:)' * expandvar   - expandvar' * sumstatevariable(1,:);
            negHz = 0;
            for jjj=1:lqp2
                znu=xn2{1}(jjj) ;
                % note that the sqrt(pi)^N+1  and 2^N factors from the integration over priors is included in the quadrature routines.
                % this makes is easier to switch between Gaussian and Uniform RV
                negHz = negHz + wn2(jjj) * (wn(:)' * log(exp(-(znu + diffsumm).^2/2/signu^2 - log(signu) -log(2*pi)/2   ) * wn(:)));
            end
    end
    % MI = H(z) - H(z|P)
    %  H(z|P)  constant ==> max MI = max H(z) = min -H(z)
    MIGaussObj = negHz;

    %%
    % Create an optimization problem using these converted optimization expressions.
    disp('create optim prob')
    convprob = optimproblem('Objective',MIGaussObj , "Constraints",stateconstraint);
    myidx = varindex(convprob )
    %%
    % View the new problem.
    %show(convprob)
    problem = prob2struct(convprob,'ObjectiveFunctionName','reducedObjective','ConstraintFunctionName','reducedConstraint');


    % compare walker solution at qp
    switch (NumberUncertain)
        case(3)
            optparams = struct('t0',[t0qp(1);0],'gammaPdfA',[alphamean(1)  ;1],'gammaPdfB',[betamean(1);1],...
                'scaleFactor',VIF_scale_fact,'T1s',[T1pmean(1),T1lmean(1)],'ExchangeTerms',[0,kplqp(1) ;0,0],...
                'TRList',TR_list,'PerfusionTerms',[kveqp(1),0],'volumeFractions',ve,...
                'fitOptions', opts);
        case(4)
            error("WIP")
        case(5)
            error("WIP")
    end
    % solve
    switch (myoptions.Algorithm)
        case('constDirect')
            powergrid = [0:1:35];
            brutesearch= zeros(size(powergrid ));
            backspaces = '';
            for iii = 1:length(powergrid(:))
                %disp(sprintf('%d %d ',pyrgrid(iii),lacgrid(iii)));
                % Print percentage progress
                percentage = iii/length(powergrid(:));
                perc_str = sprintf('completed %3.1f', percentage);
                fprintf([backspaces, perc_str]);
                backspaces = repmat(sprintf('\b'), 1, length(perc_str));

                x0.powerList = repmat(powergrid(iii),1,Ntime);
                x0.state  = evaluate(auxvariable ,x0);
                brutesearch(iii) = evaluate(MIGaussObj,x0);
            end
            save(sprintf('brutesearchNG%dNu%d%s%sSNR%02d%s.mat',NGauss,NumberUncertain,myoptions.Algorithm,ObjectiveType,modelSNR,QuadratureRule) ,'brutesearch','powergrid')
            [maxMI,idmax] = max(brutesearch(:));
            [fval,idmin] = min(brutesearch(:));
            popt.powerList = repmat(powergrid(idmin),1,Ntime);
            %pneg.FaList = repmat([pi/180*pyrgrid(idmax);pi/180*lacgrid(idmax)],1,Ntime);
            %optparams.FaList = repmat([pi/180*pyrgrid(idmin);pi/180*lacgrid(idmin)],1,Ntime);
            %refparams.FaList = repmat([pi/180*pyrgrid(idmin);pi/180*lacgrid(idmin)],1,Ntime);

            handle = figure(5)
            plot(powergrid,brutesearch)
            xlabel('power ')
            ylabel('MI ')
            title(sprintf('max %f min %f',maxMI,fval) )
            %text(pyrgrid(idmin)+1,lacgrid(idmin)+1, sprintf('opt %d %d', pyrgrid(idmin), lacgrid(idmin)));
            %text(pyrgrid(idmax)+1,lacgrid(idmax)+1, 'control');
        otherwise
            InitialGuess =  [flips(:)];
            pmin =  [flips(:)*0];
            pmax =  [flips(:)*0+35*pi/180];
            tolx=1.e-9;
            tolfun=5.e-4;
            maxiter=400;

            % truthconstraint = infeasibility(stateconstraint,x0);
            %[popt,fval,exitflag,output] = solve(convprob,x0,'Options',myoptions, 'ConstraintDerivative', 'auto-reverse', 'ObjectiveDerivative', 'auto-reverse' )
            Fx = @(x) MIGHQuadHPTofts(x, problem, myidx,Nspecies,Ntime,auxvariable);
            [designopt,fval,exitflag,output,lambda,grad,hessian] ...
                =fmincon(Fx, InitialGuess ,[],[],[],[],pmin,pmax,[],myoptions);

            handle = figure(5)
            optparams.FaList = reshape(designopt(:),size(params.FaList ));
            refparams.FaList = reshape(designopt(:),size(params.FaList ));
            popt.FaList      = reshape(designopt(:),size(params.FaList ));
    end
    % save convergence history
    set(gca,'FontSize',16)
    saveas(handle,sprintf('historyNG%dNu%d%s%sSNR%02d%s',NGauss,NumberUncertain,myoptions.Algorithm,ObjectiveType,modelSNR,QuadratureRule ),'png')
    popt.state       = evaluate(auxvariable, popt);
    toc;


    %% [t_axisopt,Mxyopt,Mzopt] = model.compile(M0.',optparams);
    %% [t_axisref,Mxyref,Mzref] = model.compile(M0.',refparams);
    %% save(sprintf('poptNG%dNu%d%s%sSNR%02d%s.mat',NGauss,NumberUncertain,myoptions.Algorithm,ObjectiveType,modelSNR,QuadratureRule) ,'fval','popt','params','Mxy','Mz','Mxyref','Mzref','signu','signuImage')
    %% handle = figure(10)
    %% plot(params.TRList,Mxyopt(1,:),'b',params.TRList,Mxyopt(2,:),'k')
    %% ylabel('MI Mxy')
    %% xlabel('sec'); legend('Pyr','Lac')
    %% set(gca,'FontSize',16)
    %% saveas(handle,sprintf('OptMxyNG%dNu%d%s%sSNR%02d%s',NGauss,NumberUncertain,myoptions.Algorithm,ObjectiveType,modelSNR,QuadratureRule),'png')
    %% handle = figure(11)
    %% plot(params.TRList,popt.FaList(1,:)*180/pi,'b',params.TRList,popt.FaList(2,:)*180/pi,'k')
    %% ylabel('MI FA (deg)')
    %% xlabel('sec'); legend('Pyr','Lac')
    %% ylim([0 40])
    %% set(gca,'FontSize',16)
    %% saveas(handle,sprintf('OptFANG%dNu%d%s%sSNR%02d%s',NGauss,NumberUncertain,myoptions.Algorithm,ObjectiveType,modelSNR,QuadratureRule),'png')
    %% handle = figure(12)
    %% plot(params.TRList,Mzopt(1,:),'b--',params.TRList,Mzopt(2,:),'k--')
    %% hold
    %% plot(params.TRList,cos(optparams.FaList(1,:))'.*popt.state(:,1, 1),'b',params.TRList,cos(optparams.FaList(2,:))'.*popt.state(:,2, 1),'k')
    %% %if(lqp > 1)
    %% %  plot(params.TRList,popt.state(:,1, 5),'b',params.TRList,popt.state(:,2, 5),'k')
    %% %  plot(params.TRList,popt.state(:,1,10),'b',params.TRList,popt.state(:,2,10),'k')
    %% %  plot(params.TRList,popt.state(:,1,15),'b',params.TRList,popt.state(:,2,15),'k')
    %% %end
    %% ylabel('MI Mz ')
    %% xlabel('sec'); legend('Pyr','Lac')
    %% set(gca,'FontSize',16)
    %% saveas(handle,sprintf('OptMzNG%dNu%d%s%sSNR%02d%s',NGauss,NumberUncertain,myoptions.Algorithm,ObjectiveType,modelSNR,QuadratureRule),'png')
end


function [MIobjfun, MIobjfun_Der]=MIGHQuadHPTofts(xopt,problem,myidx,Nspecies,Ntime,auxvariable)
x0.FaList = reshape(xopt,Nspecies,Ntime);
x0.state  = evaluate(auxvariable ,x0);
Xfull = [ x0.FaList(:); x0.state(:)];
[MIobjfun,initVals.g] = problem.objective(Xfull);
[initConst.ineq,initConst.ceq, initConst.ineqGrad,initConst.ceqGrad] = problem.nonlcon(Xfull);
objectiveGradFA    = initVals.g(myidx.FaList);
objectiveGradState = initVals.g(myidx.state);
jacobianFA    = initConst.ceqGrad(myidx.FaList,:);
jacobianState = initConst.ceqGrad(myidx.state,:);
adjointvar =-jacobianState \objectiveGradState ;
MIobjfun_Der = objectiveGradFA +  jacobianFA *   adjointvar ;
end

function tempqoi = pennesmht(k_t,w_t,rhocp_t,Md,H)
%%% Following are the inputs for the function
% rhocp_t = Tissue Density [kg/m3] * cp_t = Tissue Specific Heat [J/kg/K]
% k_t = Tissue Thermal Conductivity [W/m/K]
% w_t = Blood Perfusion Rate [1/s]
% Qm_t = Metabolic Heat Generation Rate [W/m3]
% H = Applied Magnetic Field [A/m]
% Md = Domain Magnetization Value of MNP [A/m]
% Keff = Magnetic Anistropy Constant of MNP [J/m3]
Qm_t = 0.0;
Keff = 1.0;

%%% Domain Parameters
RT = 0.005;               % [m] Radius of the Tumor
R = 5*RT;                 % [m] Radius of the computaional domain
vol_t = (4*pi*RT^3)/3;    % [m3] Tumor Volume
dr = RT/50;               % [m] Spatial discretization step size
r = (0:dr:R)';            % [m] Spatial loactions
N = length(r);            % [-] Number of spatial nodes
t_loc = r<=RT;            % [-] Tumor r indices

%%% Time discretization: Implicit method is employed
dt = 10;             % Time step [s]
t_end = 1200;       % End time [s]
t = 0:dt:t_end;     % Discrete times
TS = length(t);     % Number of time steps

%%% Ambient Conditions
htc = 10;
T_amb = 24;

%%% Blood Properties
rho_b = 1000;           % [kg/m3] Blood Density
cp_b = 3840;            % [J/kg/K] Blood Specific Heat
k_b = 0.50;             % [W/m/K] Blood Thermal Conductivity
T_b = 37;               % [°C] Blood Arterial Temperature
T_bK = T_b + 273.15;    % [K] Blood Arterial Temperature

T_initial = 37;     % [°C] Initial Condition Temperature
alpha_t = k_t/(rhocp_t); % Tissue Thermal diffusivity [m^2/s]
f = 1e5;    % [Hz] Magnetic Field Frequency f = 100 kHz

%%% Magnetic Nanoparticle Physical and Magnetic Parameters MNP Used: Fe3O4
d_mnp = 1.6e-8;                 % [m] MNP Diameter d = 16 nm
delta_l = 1e-9;                 % [m] Liquid Layer Thickness on the Hard Solid MNP delta = 1 nm
dh_mnp = d_mnp+2*delta_l;       % [m] MNP Hydrodynamic Diameter
rho_mnp = 5180;                 % [kg/m3] MNP Density

%%% MNP Dose and Magnetic Fluid Parameters
MNP_dose = 4;                   % [kg/m3] MNP Dose MNP mass used per unit vol of tumor = 4 mg/cm3; % Conversion 1 mg/cm3 = 1 kg/m3
MNP_mass = MNP_dose*vol_t;      % [kg] MNP mass injected into the tumor; Varies with MNP Dose
MNP_vol = MNP_mass/rho_mnp;     % [m3] volume of MNP in the Magnetic Fluid
MF_inj_vol = 1;                 % [ml, milliliters] Magnetic Fluid Volume injected into the Tumor
MF_vol = MF_inj_vol*1e-6;       % [m3] Magnetic Fluid Volume Injected in to the Tumor % Converion Factor for ml to m3 is 1e-6
mu_cf = 1e-3;                   % [Pa.s] Viscocity of the Carrier Fluid
MNP_vol_frac = MNP_vol/MF_vol;  %[-] Volume Fraction of MNP in MF

%%% Constants used in SAR Calculation
mu0 = pi*4e-7;              % [H/m] Permeability of Free Space
kB = 1.38e-23;              % [J/K] Blotzmann Constant
MNP_svol = pi*(d_mnp^3)/6;  % [m3] Solid Volume of Individual Mangetic Nanoparticle
MNP_hvol = pi*(dh_mnp^3)/6; % [m3] Hydrodynamic Volume of Individual Mangetic Nanoparticle
tau0 = 1e-9;                % [s] Attempt Time for Mangetic Moment relaxation
omega = 2*pi*f;             % [rad/s] Angular frequency of applied magnetic field
alpha_CF = 0.55;            % [-] Correction Factor for MNP Heating in MF and Tissue Mediums

%%% SAR Calculations
gamma = Keff*MNP_svol/(kB*T_bK);                    % [-] An intermediate parameters used in further calculations
tauB = (3*mu_cf*MNP_hvol)/(kB*T_bK);                % [s] Brownian Relaxation Time of MNP
tauN = sqrt(pi)*tau0*exp(gamma)/(2*sqrt(gamma));    % [s] Neel Relaxation Time of MNP
tauE = tauB*tauN/(tauB+tauN);                       % [s] Effective Relaxation Time of MNP
zeta = mu0*Md*H*MNP_svol/(kB*T_bK);                 % [-] Langevin parameter
X_i = mu0*(Md^2)*MNP_vol_frac*MNP_svol/(3*kB*T_bK); % [-] Initial Magnetic Susceptibility
X_0 = 3*X_i*(coth(zeta)-1/zeta)/zeta;               % [-] Equilibrium Magnetic Susceptibility
X_L = X_0*omega*tauE/(1+(omega*tauE)^2);            % [-] Loss Component of Magnetic Susceptibility
P = pi*mu0*f*X_L*H^2;                               % [W/m3] Heat Generation rate of MNPs in MF
SAR = P/(rho_mnp*MNP_vol_frac);                     % [W/kg] MNP Specific Absorption Rate

%%% Heat Source by MNP in Tumor
MNP_conc = MNP_mass/vol_t;      % [kg/m3] Concentration of MNP in Tumor, Assuming MNPs are distributed uniformly and confined within the tumor only
Q_MNP = MNP_conc*SAR;  % [W/m3] Heat Generation by MNP in Tissue
%% Implicit Finite Difference Scheme using the Cauchy Boundary Condition on the outer boundary
T = zeros(N, TS);
T(:,1) = T_initial; % Initial Condition
%%% Coefficient Matrix Three-Diagonal Generation
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
q_mnp = Q_MNP*t_loc;            % MNP heat generation [W/m^3] On for all simulation time
%%% Solver Time Iterations
for n = 2:TS
    Force = zeros(N,1); % Force vector
    Force(1:N-1)  = T(1:N-1,n-1) + (rho_b*cp_b.*w_t*dt*T_b + Qm_t*dt+ q_mnp(1:N-1)*dt)/(rhocp_t);
    Force(N) = htc*T_amb;
    T(:,n) = thomas_algorithm(Lower, Main, Upper, Force); % Spatio-Temporal Temperature
end

%%
T_sum_time = sum(T,2);      % Summation of Temperature over all Times at each r location
G = trapz(r, 4*pi*r.^2.*T_sum_time);     % Integration of Temperature over the domain in 3D
tempqoi = G;
%% Tri-Diagonal Matrix Alogorithm for Linear System of Equation Solver
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
end
