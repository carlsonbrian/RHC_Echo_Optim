% ***********************************************************************************
%             S M I T H   C A R D I O V A S C U L A R   S Y S T E M S  
%                 M O D E L   O B J E C T I V E   F U N C T I O N
% ***********************************************************************************
%
%   This function caclulates the error between the RHC data and Echo data or 
%   just the RHC data only and the Smith et al. model with adjusted parameter 
%   values. This function is called by fmincon or ga in the main script to 
%   optimize the set of adjustable parameters, p.
%
%   Model originally created on     14 November 2016
%   Model last modfied on           14 December 2018

%   Developed by        Brian Carlson
%                       Physiological Systems Dynamics Laboratory
%                       Department of Molecular and Integrative Physiology
%                       University of Michigan
%
% ***********************************************************************************
%  Start of             S M I T H   C V   M O D E L   O B J   F U N C T I O N
% ***********************************************************************************

%% **********************************************************************************
%  Optim Data for       S M I T H   C V   M O D E L   O B J   F U N C T I O N
% ***********************************************************************************

function Res = PatSpecHF_ObjFun_MinRuns(p,AllStruct_Struct)

    warning('off','all')

%%  Unpack all the passed structures
    if (numel(fieldnames(AllStruct_Struct)) == 6)
        FlagData_Struct = AllStruct_Struct.FlagData_Struct;
        PatData_Struct = AllStruct_Struct.PatData_Struct;
        RHCData_Struct = AllStruct_Struct.RHCData_Struct;
        EchoData_Struct = AllStruct_Struct.EchoData_Struct;
        CVParam_Struct = AllStruct_Struct.CVParam_Struct;
        SimOptParam_Struct = AllStruct_Struct.SimOptParam_Struct;
    else 
        FlagData_Struct = AllStruct_Struct.FlagData_Struct;
        PatData_Struct = AllStruct_Struct.PatData_Struct;
        RHCData_Struct = AllStruct_Struct.RHCData_Struct;
        CVParam_Struct = AllStruct_Struct.CVParam_Struct;
        SimOptParam_Struct = AllStruct_Struct.SimOptParam_Struct;
    end
        
%%  Unpack the RHC and Echo data or RHC data alone
    
    RHCEcho_Flag = FlagData_Struct.RHCEcho_Flag;
    if (RHCEcho_Flag == 1)
        % Unpack patient data
        BW_Data = PatData_Struct.BW;                % Body weight (kg)
        Hgt_Data = PatData_Struct.Hgt;              % Height (cm)
        Gender_Data = PatData_Struct.Gender;        % Gender (M or F)
        % Unpack RHC data
        P_RVsyst_Data = RHCData_Struct.P_RVsyst;    % Syst RV pressure (mmHg)
        P_RVdiast_Data = RHCData_Struct.P_RVdiast;  % Diast RV pressure (mmHg)
        P_PAsyst_Data = RHCData_Struct.P_PAsyst;    % Syst pulm art press (mmHg)
        P_PAdiast_Data = RHCData_Struct.P_PAdiast;  % Diast pulm art press (mmHg)
        P_PCWave_Data = RHCData_Struct.P_PCWave;    % Ave pulm wedge press (mmHg)
        P_AOsyst_Data = RHCData_Struct.P_AOsyst;    % Syst aortic press (mmHg)
        P_AOdiast_Data = RHCData_Struct.P_AOdiast;  % Diast aortic press (mmHg)
        HR_RHC_Data = RHCData_Struct.HR_RHC;        % RHC heart rate (beats/min)
        CO_Fick_Data = RHCData_Struct.CO_Fick;      % Cardiac output Fick (L/min)
        CO_Thermo_Data = RHCData_Struct.CO_Thermo;  % Crd out thrmdil (L/min)
        % Unpack Echo data
        V_LVsyst_Data = EchoData_Struct.V_LVsyst;   % Systolic LV volume (mL)   
        V_LVdiast_Data = EchoData_Struct.V_LVdiast; % Diastolic LV volume (mL)
        HR_Echo_Data   = EchoData_Struct.HR_Echo;   % Echo heart rate (beats/min)
        CO_EchoD_Data = EchoData_Struct.CO_EchoD;   % Crdc outpt Echo-Dop (L/min)
        
    else
        % Unpack patient data
        BW_Data = PatData_Struct.BW;                % Body weight (kg)
        Hgt_Data = PatData_Struct.Hgt;              % Height (cm)
        Gender_Data = PatData_Struct.Gender;        % Gender (M or F)
        % Unpack RHC data
        P_RVsyst_Data = RHCData_Struct.P_RVsyst;    % Syst RV pressure (mmHg)
        P_RVdiast_Data = RHCData_Struct.P_RVdiast;  % Diast RV pressure (mmHg)
        P_PAsyst_Data = RHCData_Struct.P_PAsyst;    % Syst pulm art press (mmHg)
        P_PAdiast_Data = RHCData_Struct.P_PAdiast;  % Diast pulm art press (mmHg)
        P_PCWave_Data = RHCData_Struct.P_PCWave;    % Ave pulm wedge press (mmHg)
        P_AOsyst_Data = RHCData_Struct.P_AOsyst;    % Syst aortic press (mmHg)
        P_AOdiast_Data = RHCData_Struct.P_AOdiast;  % Diast aortic press (mmHg)
        HR_RHC_Data = RHCData_Struct.HR_RHC;        % Ave heart rate (beats/min)
        CO_Fick_Data = RHCData_Struct.CO_Fick;      % Cardiac output Fick (L/min)
        CO_Thermo_Data = RHCData_Struct.CO_Thermo;  % Crd out thrmdil (L/min)
    end
    
    
%%  Unpack the simulation parameters 
    AdjParam_Strngs = SimOptParam_Struct.AdjParam_Strngs;
    NumBeats_SS = SimOptParam_Struct.NumBeats_SS;
    NumBeats_ResPlot = SimOptParam_Struct.NumBeats_ResPlot;
    
    
%%  Unpack the adjustable parameters and write into CVParam structure

    Num_AdjParams = size(p,2);
    CVParam_FieldNames = fieldnames(CVParam_Struct);
    Num_AllParams = size(CVParam_FieldNames,1);
    for i = 1:Num_AdjParams
        for j = 1:Num_AllParams
            if (strcmp(CVParam_FieldNames{j},AdjParam_Strngs{i}))
                if (strcmp(AdjParam_Strngs{i},'P_th'))
                    CVParam_Struct.(AdjParam_Strngs{i}) = p(i);
                else
                    CVParam_Struct.(AdjParam_Strngs{i}) = exp(p(i));
                end
            end
        end
    end
    

%%  Calculate circulating blood volume for all simulations

    SVFact = CVParam_Struct.SVFact;
    % Calculate total blood volume based on height, weight and gender.
    %  This expression is from Nadler et al. Surgery 51:224,1962.
    if (Gender_Data == 'M')
        TotBV = ((0.3669 * (Hgt_Data/100)^3) + (0.03219 * BW_Data) + 0.6041) * 1000;
    else
        TotBV = ((0.3561 * (Hgt_Data/100)^3) + (0.03308 * BW_Data) + 0.1833) * 1000;
    end
    % The original Smith model only circulated a portion of the blood so aortic
    %  pressure dynamics are not lumped into a general arterial systemic 
    %  compartment. Assuming they were simulating a typical 5000 mL total blood
    %  volume they included only 1500 mL (or 30%) in the circulating volume
    %  therefore we will multiply our calculated TotBV value by 0.3 to yield 
    %  circulating blood volume. To account for extra recruited volume in heart
    %  disease the 30% circulating blood volume can be altered by changing SVFact
    CircBV = SVFact * 0.30 * TotBV;
    
    % Setting state variable initial conditions 
    %  Note that initial volume division is scaled as a fraction
    %  of circulating blood volume calculated earlier and the 
    %  initial guess of flow is total blood volume in one minute
    V_lv0 = (94.6812/1500) * CircBV;
    V_rv0 = (90.7302/1500) * CircBV;
    V_pa0 = (43.0123/1500) * CircBV;
    V_pu0 = (808.458/1500) * CircBV;
    V_ao0 = (133.338/1500) * CircBV;
    V_vc0 = (329.780/1500) * CircBV;
    Q_mt0 = TotBV/60;
    Q_av0 = TotBV/60;
    Q_tc0 = TotBV/60;
    Q_pv0 = TotBV/60;
    
    
%%  Execute simulations for RHC only or RHC/Echo data

    if (RHCEcho_Flag == 1)                          % If 1 both RHC and Echo data
        
        % Use a try catch in case the parameter values cause a crash 
        try
            
            % Set initial conditions on explicit state
            %  variables for both RHC and Echo simulations
            X0(1) = V_lv0;
            X0(2) = V_rv0;
            X0(3) = V_pa0;
            X0(4) = V_pu0;
            X0(5) = V_ao0;
            X0(6) = V_vc0;
            X0(7) = Q_mt0;
            X0(8) = Q_av0;
            X0(9) = Q_tc0;
            X0(10) = Q_pv0;
        
            % FIRST PERFORM THE RHC SIMULATION
            % Build a structure with the RHC driver parameters
            HR = HR_RHC_Data;                       % RHC heart rate (beats/min)
            period = 60/HR_RHC_Data;                % Period of heart beat (s)
            B = HR_RHC_Data;                        % Elastance fctn param (1/s^2)
            C = period/2;                           % Elastance fctn param (s)
            DriverP_Values = {HR period B C};
            DriverP_Fields = {'HR' 'period' 'B' 'C'};
            DriverP_Struct = cell2struct(DriverP_Values, ...
                DriverP_Fields,2);
            % Calculating timespans to reach steady state
            TSpan_SS = [0 NumBeats_SS * period];
            % Now set the starting time and ventricular
            %  volumes to start the simulations
            time = 0;
            V_lv0 = X0(1);
            V_rv0 = X0(2);
            fsolveOpts = optimoptions('fsolve');
            fsolveOpts.Display = 'off';
            % Calculating the initial condition on the volume due to the septum 
            %  deflection which is an implicit function requiring the use of fsolve
            SeptZF_Hndl = @(V_spt) SeptZF(V_spt,V_lv0,V_rv0, ...
                time,DriverP_Struct,CVParam_Struct);
            [V_spt0,~,~,~] = fsolve(SeptZF_Hndl,-15,fsolveOpts);
            X0(11) = V_spt0;
            % Build mass matrix for DAE solver
            M = eye(11);                            % Put identity on diagonal
            M(11,11) = 0;                           % Set last express as a ZFun
            % Set ODE/DAE options and time span
            ODE_Opts = odeset('Mass',M); 
            % Solve over the steady state time span with ode15s
            [~,X_RHC_Out_SS] = ...
                ode15s(@dXdT_Smith,TSpan_SS,X0, ...
                ODE_Opts,DriverP_Struct,CVParam_Struct);
            % Now solve over the simulation time span
            [T_RHC_Out,X_RHC_Out] = ...
                ode15s(@dXdT_Smith,TSpan_Sim,X_RHC_Out_SS(end,:), ...
                ODE_Opts,DriverP_Struct,CVParam_Struct);
            %  Now run Smith et al. model to get intermediate pressures
            Num_TOut_RHCSim = size(T_RHC_Out,1);    % Number of time points
            P_RV_RHCSim = zeros(Num_TOut_RHCSim,1); % Preallocate matrices
            P_AO_RHCSim = zeros(Num_TOut_RHCSim,1);
            P_PA_RHCSim = zeros(Num_TOut_RHCSim,1);
            P_PU_RHCSim = zeros(Num_TOut_RHCSim,1);
            % Use output state variable vector to get pressures
            for i = 1:Num_TOut_RHCSim
                VarOut = dXdT_Smith(T_RHC_Out(i), ...
                    X_RHC_Out(i,:),DriverP_Struct,CVParam_Struct,1);
                P_RV_RHCSim(i) = VarOut(2);
                P_AO_RHCSim(i) = VarOut(3);
                P_PA_RHCSim(i) = VarOut(5);
                P_PU_RHCSim(i) = VarOut(6);
            end
            
            % GET SIMULATION VALUES TO COMPARE TO RHC DATA
            % Assigning initial low vlaues on systolic and 
            %  high values on diastolic pressures
            P_RVsyst_RHCSim = 0;
            P_RVdiast_RHCSim = 200;
            P_PAsyst_RHCSim = 0;
            P_PAdiast_RHCSim = 150;
            P_AOsyst_RHCSim = 0;
            P_AOdiast_RHCSim = 250;
            V_LVsyst_RHCSim = 150;
            V_LVdiast_RHCSim = 0;
            P_PCWave_RHCSum = 0;
            % Specifying the number of beats at the end of the RHC simulation
            %  that we want to calculate the residual from and then finding 
            %  the portion of the simulation to extract values from
            TimeStart_RHCResCalc = ...
                (NumBeats_Sim - NumBeats_PlotRes) * period;
            tIndStart_RHCResCalc = ...
                find((T_RHC_Out >= TimeStart_RHCResCalc),1,'first');
            % Looking for systolic and diastolic values in simulation results
            for i = tIndStart_RHCResCalc:Num_TOut_RHCSim
                P_RVsyst_RHCSim = max(P_RVsyst_RHCSim,P_RV_RHCSim(i));
                P_RVdiast_RHCSim = min(P_RVdiast_RHCSim,P_RV_RHCSim(i));
                P_PAsyst_RHCSim = max(P_PAsyst_RHCSim,P_PA_RHCSim(i));
                P_PAdiast_RHCSim = min(P_PAdiast_RHCSim,P_PA_RHCSim(i));
                P_AOsyst_RHCSim = max(P_AOsyst_RHCSim,P_AO_RHCSim(i));
                P_AOdiast_RHCSim = min(P_AOdiast_RHCSim,P_AO_RHCSim(i));
                P_PCWave_RHCSum = P_PCWave_RHCSum + P_PU_RHCSim(i);
                V_LVsyst_RHCSim = min(V_LVsyst_RHCSim,X_RHC_Out(i,1));
                V_LVdiast_RHCSim = max(V_LVdiast_RHCSim,X_RHC_Out(i,1));
            end
            P_PCWave_RHCSim = P_PCWave_RHCSum / Num_TOut_RHCSim;
            CO_RHCSim = ((V_LVdiast_RHCSim - V_LVsyst_RHCSim) * HR_RHC_Data) / 1000;
            
            % NEXT PERFORM THE ECHO SIMULATION
            % Build a structure with the Echo driver parameters
            HR = HR_Echo_Data;                      % RHC heart rate (beats/min)
            period = 60/HR_Echo_Data;               % Period of heart beat (s)
            B = HR_Echo_Data;                       % Elastance fctn param (1/s^2)
            C = period/2;                           % Elastance fctn param (s)
            DriverP_Values = {HR period B C};
            DriverP_Fields = {'HR' 'period' 'B' 'C'};
            DriverP_Struct = cell2struct(DriverP_Values, ...
                DriverP_Fields,2);
            % Now set the starting time and ventricular
            %  volumes to start the simulations
            time = 0;
            V_lv0 = X0(1);
            V_rv0 = X0(2);
            fsolveOpts = optimoptions('fsolve');
            fsolveOpts.Display = 'off';
            % Calculating the initial condition on the volume due to the septum 
            %  deflection which is an implicit function requiring the use of fsolve
            SeptZF_Hndl = @(V_spt) SeptZF(V_spt,V_lv0,V_rv0, ...
                time,DriverP_Struct,CVParam_Struct);
            [V_spt0,~,~,~] = fsolve(SeptZF_Hndl,-15,fsolveOpts);
            X0(11) = V_spt0;
            % Build mass matrix for DAE solver
            M = eye(11);                            % Put identity on diagonal
            M(11,11) = 0;                           % Set last express as a ZFun
            % Set ODE/DAE options and time span
            ODE_Opts = odeset('Mass',M); 
            % Solve over the steady state time span with ode15s
            [~,X_Echo_Out_SS] = ...
                ode15s(@dXdT_Smith,TSpan_SS,X0, ...
                ODE_Opts,DriverP_Struct,CVParam_Struct);
            % Now solve over the simulation time span
            [T_Echo_Out,X_Echo_Out] = ...
                ode15s(@dXdT_Smith,TSpan_Sim,X_Echo_Out_SS(end,:), ...
                ODE_Opts,DriverP_Struct,CVParam_Struct);
            
            % GET SIMULATION VALUES TO COMPARE TO ECHO DATA
            % Assigning initial low vlaues on systolic and 
            %  high values on diastolic pressures
            Num_TOut_EchoSim = size(T_Echo_Out,1);  % Number of time points
            V_LVsyst_EchoSim = 150;
            V_LVdiast_EchoSim = 0;
            % Specifying the number of beats at the end of the Echo simulation
            %  that we want to calculate the residual from and then finding
            %  the portion of the simulation to extract values from
            TimeStart_EchoResCalc = ...
                (NumBeats_Sim - NumBeats_PlotRes) * period;
            tIndStart_EchoResCalc = ...
                find((T_RHC_Out >= TimeStart_EchoResCalc),1,'first');
            % Looking for systolic and diastolic values in simulation results
            for i = tIndStart_EchoResCalc:Num_TOut_EchoSim
                V_LVsyst_EchoSim = min(V_LVsyst_EchoSim,X_Echo_Out(i,1));
                V_LVdiast_EchoSim = max(V_LVdiast_EchoSim,X_Echo_Out(i,1));
            end
            % Calculate Echo CO
            CO_EchoSim = ((V_LVdiast_EchoSim - ...
                V_LVsyst_EchoSim) * HR_Echo_Data) / 1000;
            
            % NOW CALCULATING THE RESIDUAL FOR BOTH RHC AND ECHO DATA
            %  Calculate the difference between data and simulation
            % Right ventricle pressure normalized residual
            P_RVsyst_Res = abs(P_RVsyst_RHCSim - P_RVsyst_Data) / P_AOsyst_Data;
            P_RVdiast_Res = abs(P_RVdiast_RHCSim - P_RVdiast_Data) / P_AOsyst_Data;
            % Pulmonary artery pressure normalized residual
            P_PAsyst_Res = abs(P_PAsyst_RHCSim - P_PAsyst_Data) / P_AOsyst_Data;
            P_PAdiast_Res = abs(P_PAdiast_RHCSim - P_PAdiast_Data) / P_AOsyst_Data;
            % Aortic pressure normalized residual
            P_AOsyst_Res = abs(P_AOsyst_RHCSim - P_AOsyst_Data) / P_AOsyst_Data;
            P_AOdiast_Res = abs(P_AOdiast_RHCSim - P_AOdiast_Data) / P_AOsyst_Data;
            % Pulmonary wedge pressure (in pulmonary vein) normalized residual
            P_PCWave_Res = abs(P_PCWave_RHCSim - P_PCWave_Data) / P_AOsyst_Data;
            % Cardiac output normalized residual for both simulations
            CO_EchoSVHR = ((V_LVdiast_Data - V_LVsyst_Data) * HR_Echo_Data) / 1000;
            CO_RHCEcho_Max = max([CO_Thermo_Data CO_EchoD_Data CO_EchoSVHR]);
            CO_RHCRes = abs(CO_RHCSim - CO_Thermo_Data) / CO_RHCEcho_Max;
            CO_EchoSVHRRes = abs(CO_EchoSim - CO_EchoSVHR) / CO_RHCEcho_Max;
            CO_EchoDRes = abs(CO_EchoSim - CO_EchoD_Data) / CO_RHCEcho_Max;
            CO_EchoRes = (CO_EchoSVHRRes + CO_EchoDRes) / 2;
            % Left ventricular volume normalized residuals
            V_LVsyst_Res = abs(V_LVsyst_EchoSim - ...
                V_LVsyst_Data) / V_LVdiast_Data;
            V_LVdiast_Res = abs(V_LVdiast_EchoSim - ...
                V_LVdiast_Data) / V_LVdiast_Data;
    
            Res = (P_RVsyst_Res + P_RVdiast_Res + P_PAsyst_Res + ...
                P_PAdiast_Res + P_AOsyst_Res + P_AOdiast_Res + ...
                P_PCWave_Res + CO_RHCRes + CO_EchoRes + ...
                V_LVsyst_Res + V_LVdiast_Res) / 11;
            
        catch
            
            Res = 10;
            
        end
        
    else                                            % Only RHC data
         
        try
            
            % Set initial conditions on explicit state
            %  variables for RHC simulation
            X0(1) = V_lv0;
            X0(2) = V_rv0;
            X0(3) = V_pa0;
            X0(4) = V_pu0;
            X0(5) = V_ao0;
            X0(6) = V_vc0;
            X0(7) = Q_mt0;
            X0(8) = Q_av0;
            X0(9) = Q_tc0;
            X0(10) = Q_pv0;
        
            % PERFORM THE RHC SIMULATION
            % Build a structure with the RHC driver parameters
            HR = HR_RHC_Data;                       % RHC heart rate (beats/min)
            period = 60/HR_RHC_Data;                % Period of heart beat (s)
            B = HR_RHC_Data;                        % Elastance fctn param (1/s^2)
            C = period/2;                           % Elastance fctn param (s)
            DriverP_Values = {HR period B C};
            DriverP_Fields = {'HR' 'period' 'B' 'C'};
            DriverP_Struct = cell2struct(DriverP_Values, ...
                DriverP_Fields,2);
            % Calculating timespans to reach steady state and then for simulation
            TSpan_SS = [0 NumBeats_SS * period];
            TSpan_Sim = [0 NumBeats_Sim * period];
            % Now set the starting time and ventricular
            %  volumes to start the simulations
            time = 0;
            V_lv0 = X0(1);
            V_rv0 = X0(2);
            fsolveOpts = optimoptions('fsolve');
            fsolveOpts.Display = 'off';
            % Calculating the initial condition on the volume due to the septum 
            %  deflection which is an implicit function requiring the use of fsolve
            SeptZF_Hndl = @(V_spt) SeptZF(V_spt,V_lv0,V_rv0, ...
                time,DriverP_Struct,CVParam_Struct);
            [V_spt0,~,~,~] = fsolve(SeptZF_Hndl,-15,fsolveOpts);
            X0(11) = V_spt0;
            % Build mass matrix for DAE solver
            M = eye(11);                            % Put identity on diagonal
            M(11,11) = 0;                           % Set last express as a ZFun
            % Set ODE/DAE options and time span
            ODE_Opts = odeset('Mass',M); 
            % Solve over the steady state time span with ode15s
            [~,X_RHC_Out_SS] = ...
                ode15s(@dXdT_Smith,TSpan_SS,X0, ...
                ODE_Opts,DriverP_Struct,CVParam_Struct);
            % Now solve over the simulation time span
            [T_RHC_Out,X_RHC_Out] = ...
                ode15s(@dXdT_Smith,TSpan_Sim,X_RHC_Out_SS(end,:), ...
                ODE_Opts,DriverP_Struct,CVParam_Struct);
            %  Now run Smith et al. model to get intermediate pressures
            Num_TOut_RHCSim = size(T_RHC_Out,1);    % Number of time points
            P_RV_RHCSim = zeros(Num_TOut_RHCSim,1); % Preallocate matrices
            P_AO_RHCSim = zeros(Num_TOut_RHCSim,1);
            P_PA_RHCSim = zeros(Num_TOut_RHCSim,1);
            P_PU_RHCSim = zeros(Num_TOut_RHCSim,1);
            % Use output state variable vector to get pressures
            for i = 1:Num_TOut_RHCSim
                VarOut = dXdT_Smith(T_RHC_Out(i), ...
                    X_RHC_Out(i,:),DriverP_Struct,CVParam_Struct,1);
                P_RV_RHCSim(i) = VarOut(2);
                P_AO_RHCSim(i) = VarOut(3);
                P_PA_RHCSim(i) = VarOut(5);
                P_PU_RHCSim(i) = VarOut(6);
            end
            
            % GET SIMULATION VALUES TO COMPARE TO DATA
            % Assigning initial low vlaues on systolic and 
            %  high values on diastolic pressures
            P_RVsyst_RHCSim = 0;
            P_RVdiast_RHCSim = 200;
            P_PAsyst_RHCSim = 0;
            P_PAdiast_RHCSim = 150;
            P_AOsyst_RHCSim = 0;
            P_AOdiast_RHCSim = 250;
            V_LVsyst_RHCSim = 150;
            V_LVdiast_RHCSim = 0;
            P_PCWave_RHCSum = 0;
            % Specifying the number of beats at the end of the RHC simulation
            %  that we want to calculate the residual from and then finding 
            %  the portion of the simulation to extract values from
            TimeStart_RHCResCalc = ...
                (NumBeats_Sim - NumBeats_PlotRes) * period;
            tIndStart_RHCResCalc = ...
                find((T_RHC_Out >= TimeStart_RHCResCalc),1,'first');
            % Looking for systolic and diastolic values in simulation results
            for i = tIndStart_RHCResCalc:Num_TOut_RHCSim
                P_RVsyst_RHCSim = max(P_RVsyst_RHCSim,P_RV_RHCSim(i));
                P_RVdiast_RHCSim = min(P_RVdiast_RHCSim,P_RV_RHCSim(i));
                P_PAsyst_RHCSim = max(P_PAsyst_RHCSim,P_PA_RHCSim(i));
                P_PAdiast_RHCSim = min(P_PAdiast_RHCSim,P_PA_RHCSim(i));
                P_AOsyst_RHCSim = max(P_AOsyst_RHCSim,P_AO_RHCSim(i));
                P_AOdiast_RHCSim = min(P_AOdiast_RHCSim,P_AO_RHCSim(i));
                P_PCWave_RHCSum = P_PCWave_RHCSum + P_PU_RHCSim(i);
                V_LVsyst_RHCSim = min(V_LVsyst_RHCSim,X_RHC_Out(i,1));
                V_LVdiast_RHCSim = max(V_LVdiast_RHCSim,X_RHC_Out(i,1));
            end
            P_PCWave_RHCSim = P_PCWave_RHCSum / Num_TOut_RHCSim;
            CO_RHCSim = ((V_LVdiast_RHCSim - V_LVsyst_RHCSim) * HR_RHC_Data) / 1000;
            
            % NOW CALCULATING THE RESIDUAL FOR RHC DATA
            %  Calculate the difference between data and simulation
            % Right ventricle pressure normalized residual
            P_RVsyst_Res = abs(P_RVsyst_RHCSim - P_RVsyst_Data) / P_AOsyst_Data;
            P_RVdiast_Res = abs(P_RVdiast_RHCSim - P_RVdiast_Data) / P_AOsyst_Data;
            % Pulmonary artery pressure normalized residual
            P_PAsyst_Res = abs(P_PAsyst_RHCSim - P_PAsyst_Data) / P_AOsyst_Data;
            P_PAdiast_Res = abs(P_PAdiast_RHCSim - P_PAdiast_Data) / P_AOsyst_Data;
            % Aortic pressure normalized residual
            P_AOsyst_Res = abs(P_AOsyst_RHCSim - P_AOsyst_Data) / P_AOsyst_Data;
            P_AOdiast_Res = abs(P_AOdiast_RHCSim - P_AOdiast_Data) / P_AOsyst_Data;
            % Pulmonary wedge pressure (in pulmonary vein) normalized residual
            P_PCWave_Res = abs(P_PCWave_RHCSim - P_PCWave_Data) / P_AOsyst_Data;
            % Cardiac output normalized residual
            CO_Res = abs(CO_RHCSim - CO_Thermo_Data) / CO_Thermo_Data;
    
            Res = (P_RVsyst_Res + P_RVdiast_Res + P_PAsyst_Res + ...
                P_PAdiast_Res + P_AOsyst_Res + P_AOdiast_Res + ...
                P_PCWave_Res + CO_Res) / 8;
            
        catch
            
            Res = 10;
            
        end
        
    end
         
%     %  Plot out intermediate fits  
%         if (IGTune_Flag == 0 && Parallel_Flag == 0)
%             
%             % Right ventricular pressure update
%             set(PRV_Line,'XData',T_Out,'YData',P_RVSim)
%             % Aortic pressure update
%             set(PAO_Line,'XData',T_Out,'YData',P_AOSim)
%             % Pulmonary artery pressure update
%             set(PPA_Line,'XData',T_Out,'YData',P_PASim)
%             % Pulmonary capillary wedge pressure update
%             set(PPU_Line,'XData',T_Out,'YData',P_PUSim)
%             % Simulated cardiac output and residual update
%             set(COSim_Text,'String',['CO Sim  = ' num2str(CO_Sim)])
%             set(Res_Text,'String',['Res = ', num2str(Res)])
%             % Left and right Ventricular volume update
%             set(VLV_Line,'XData',T_Out,'YData',X_Out(:,1))
%             set(VRV_Line,'XData',T_Out,'YData',X_Out(:,2))
%             drawnow
%             pause(0.15)
%             
%         end
    
    
end
