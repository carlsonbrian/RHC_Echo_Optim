% ***********************************************************************************
%             S M I T H   C A R D I O V A S C U L A R   S Y S T E M S  
%                 M O D E L   O B J E C T I V E   F U N C T I O N
% ***********************************************************************************
%
%   This function caclulates the error between the RHC data and Echo data or 
%   just the RHC data only and the reduced version of the Smith et al. model with
%   adjusted parameter values. This function is called by fmincon or ga in the 
%   main script to optimize the set of adjustable parameters, p.
%
%   Model originally created on     14 November 2016
%   Model last modfied on           16     July 2019

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

function Res = PatSpecHFRed4_ObjFun(p,AllStruct_Struct)

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
        Sex_Data = PatData_Struct.Sex;              % Sex (M or F)
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
        ID_LVsyst_Data = EchoData_Struct.ID_LVsyst;   % Systolic LV volume (mL)   
        ID_LVdiast_Data = EchoData_Struct.ID_LVdiast; % Diastolic LV volume (mL)
        HR_Echo_Data   = EchoData_Struct.HR_Echo;   % Echo heart rate (beats/min)
        CO_EchoD_Data = EchoData_Struct.CO_EchoD;   % Crdc outpt Echo-Dop (L/min)
        V_LVsyst_Data = EchoData_Struct.V_LVsyst;   % Systolic LV volume (mL)   
        V_LVdiast_Data = EchoData_Struct.V_LVdiast; % Diastolic LV volume (mL)
        
    else
        % Unpack patient data
        BW_Data = PatData_Struct.BW;                % Body weight (kg)
        Hgt_Data = PatData_Struct.Hgt;              % Height (cm)
        Sex_Data = PatData_Struct.Sex;              % Sex (M or F)
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
    % Calculate total blood volume based on height, weight and sex.
    %  This expression is from Nadler et al. Surgery 51:224,1962.
    if (Sex_Data == 'M')
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
        
            % FIRST PERFORM THE RHC SIMULATION
            % Build a structure with the RHC driver parameters
            HR = HR_RHC_Data;                       % RHC heart rate (beats/min)
            period_RHC = 60/HR_RHC_Data;            % Period of heart beat (s)
            B_RHC = HR_RHC_Data;                    % Elastance fctn param (1/s^2)
            C_RHC = period_RHC/2;                   % Elastance fctn param (s)
            DriverP_Values = {HR period_RHC B_RHC C_RHC};
            DriverP_Fields = {'HR' 'period' 'B' 'C'};
            DriverP_Struct = cell2struct(DriverP_Values, ...
                DriverP_Fields,2);
            % Calculating timespans to reach steady state
            TSpan_SS = [0 NumBeats_SS * period_RHC];
            % Solve over the steady state time span with ode15s
            [T_RHC_Out,X_RHC_Out] = ...
                ode15s(@dXdT_SmithRed4,TSpan_SS,X0, ...
                [],DriverP_Struct,CVParam_Struct);
           
            %  Now run Smith et al. model to get intermediate pressures
            Num_TOut_RHC = size(T_RHC_Out,1);    % Number of time points
            P_RV_RHC = zeros(Num_TOut_RHC,1);    % Preallocate matrices
            P_AO_RHC = zeros(Num_TOut_RHC,1);
            P_PA_RHC = zeros(Num_TOut_RHC,1);
            P_PU_RHC = zeros(Num_TOut_RHC,1);
            % Use output state variable vector to get pressures
            for i = 1:Num_TOut_RHC
                VarOut = dXdT_SmithRed4(T_RHC_Out(i), ...
                    X_RHC_Out(i,:),DriverP_Struct,CVParam_Struct,1);
                P_RV_RHC(i) = VarOut(2);
                P_AO_RHC(i) = VarOut(3);
                P_PA_RHC(i) = VarOut(5);
                P_PU_RHC(i) = VarOut(6);
            end
            
            % NEXT PERFORM THE ECHO SIMULATION
            % Build a structure with the Echo driver parameters
            HR = HR_Echo_Data;                      % RHC heart rate (beats/min)
            period_Echo = 60/HR_Echo_Data;          % Period of heart beat (s)
            B_Echo = HR_Echo_Data;                  % Elastance fctn param (1/s^2)
            C_Echo = period_Echo/2;                 % Elastance fctn param (s)
            DriverP_Values = {HR period_Echo B_Echo C_Echo};
            DriverP_Fields = {'HR' 'period' 'B' 'C'};
            DriverP_Struct = cell2struct(DriverP_Values, ...
                DriverP_Fields,2);
            % Calculating timespans to reach steady state
            TSpan_SS = [0 NumBeats_SS * period_Echo];
            
            % Solve over the steady state time span with ode15s
            [T_Echo_Out,X_Echo_Out] = ...
                ode15s(@dXdT_SmithRed4,TSpan_SS,X0, ...
                [],DriverP_Struct,CVParam_Struct);
            
            % NOW CHECK THAT THE ODES WERE INTEGRATED OVER THE FULL TIME SPAN
            T_SSRHC = NumBeats_SS * (60/HR_RHC_Data);
            T_SSEcho = NumBeats_SS * (60/HR_Echo_Data);
            if ((T_RHC_Out(end) == T_SSRHC) && ...
                    (T_Echo_Out(end) == T_SSEcho))
           
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
                    (NumBeats_SS - NumBeats_ResPlot) * period_RHC;
                tIndStart_RHCResCalc = ...
                    find((T_RHC_Out >= TimeStart_RHCResCalc),1,'first');
                % Looking for systolic and diastolic values in simulation results
                
                for i = tIndStart_RHCResCalc:Num_TOut_RHC
                    P_RVsyst_RHCSim = max(P_RVsyst_RHCSim,P_RV_RHC(i));
                    P_RVdiast_RHCSim = min(P_RVdiast_RHCSim,P_RV_RHC(i));
                    P_PAsyst_RHCSim = max(P_PAsyst_RHCSim,P_PA_RHC(i));
                    P_PAdiast_RHCSim = min(P_PAdiast_RHCSim,P_PA_RHC(i));
                    P_AOsyst_RHCSim = max(P_AOsyst_RHCSim,P_AO_RHC(i));
                    P_AOdiast_RHCSim = min(P_AOdiast_RHCSim,P_AO_RHC(i));
                    P_PCWave_RHCSum = P_PCWave_RHCSum + P_PU_RHC(i);
                    V_LVsyst_RHCSim = min(V_LVsyst_RHCSim,X_RHC_Out(i,1));
                    V_LVdiast_RHCSim = max(V_LVdiast_RHCSim,X_RHC_Out(i,1));
                end
                P_PCWave_RHCSim = P_PCWave_RHCSum / ...
                    (Num_TOut_RHC - tIndStart_RHCResCalc + 1);
                CO_RHCSim = ...
                    ((V_LVdiast_RHCSim - V_LVsyst_RHCSim) * ...
                    HR_RHC_Data) / 1000;
            
                % GET SIMULATION VALUES TO COMPARE TO ECHO DATA
                % Assigning initial low vlaues on systolic and 
                %  high values on diastolic pressures
                Num_TOut_Echo = size(T_Echo_Out,1);  % Number of time points
                V_LVsyst_EchoSim = 150;
                V_LVdiast_EchoSim = 0;
                % Specifying the number of beats at the end of the Echo simulation
                %  that we want to calculate the residual from and then finding
                %  the portion of the simulation to extract values from
                TimeStart_EchoResCalc = ...
                    (NumBeats_SS - NumBeats_ResPlot) * period_Echo;
                tIndStart_EchoResCalc = ...
                    find((T_Echo_Out >= TimeStart_EchoResCalc),1,'first');
                % Looking for systolic and diastolic values in simulation results
                for i = tIndStart_EchoResCalc:Num_TOut_Echo
                    V_LVsyst_EchoSim = min(V_LVsyst_EchoSim,X_Echo_Out(i,1));
                    V_LVdiast_EchoSim = max(V_LVdiast_EchoSim,X_Echo_Out(i,1));
                end
                % Calculate Echo CO
                CO_EchoSim = ((V_LVdiast_EchoSim - ...
                    V_LVsyst_EchoSim) * HR_Echo_Data) / 1000;
            
                % NOW CALCULATING THE RESIDUAL FOR BOTH RHC AND ECHO DATA
                %  Calculate the difference between data and simulation
                % Right ventricle pressure normalized residual
                P_RVsyst_Res = abs(P_RVsyst_RHCSim - ...
                    P_RVsyst_Data) / P_AOsyst_Data;
                P_RVdiast_Res = abs(P_RVdiast_RHCSim - ...
                    P_RVdiast_Data) / P_AOsyst_Data;
                % Pulmonary artery pressure normalized residual
                P_PAsyst_Res = abs(P_PAsyst_RHCSim - ...
                    P_PAsyst_Data) / P_AOsyst_Data;
                P_PAdiast_Res = abs(P_PAdiast_RHCSim - ...
                    P_PAdiast_Data) / P_AOsyst_Data;
                % Aortic pressure normalized residual
                P_AOsyst_Res = abs(P_AOsyst_RHCSim - ...
                    P_AOsyst_Data) / P_AOsyst_Data;
                P_AOdiast_Res = abs(P_AOdiast_RHCSim - ...
                    P_AOdiast_Data) / P_AOsyst_Data;
                % Pulmonary wedge pressure (in pulmonary vein) normalized residual
%                 P_PCWave_RHCSim
%                 P_PCWave_Data
                P_PCWave_Res = abs(P_PCWave_RHCSim - ...
                    P_PCWave_Data) / P_AOsyst_Data;
                % Cardiac output normalized residual for both simulations
                % Find normalizing CO measure
                if (V_LVsyst_Data == -1)
                    V_LVsyst_Teich = ...                    % Calculate the 
                        ((ID_LVsyst_Data/10)^3) / ...       %  LV systolic volume
                        ((6/pi) * ((0.075 * ...             %  using the Teichholz
                        (ID_LVsyst_Data/10)) + 0.18));      %  formula
                    V_LVdiast_Teich = ...                   % Calculate the 
                        ((ID_LVdiast_Data/10)^3) / ...      %  LV diastolic volume     
                        ((6/pi) * ((0.075 * ...             %  using the Teichholz    
                        (ID_LVdiast_Data/10)) + 0.18));     %  formula
                    CO_EchoTeich = ((V_LVdiast_Teich - ...  % Calculate SV * HR
                        V_LVsyst_Teich) * ...               %  estimate of cardiac
                        HR_Echo_Data) / 1000;               %  output for Teich vols
                    CO_EchoSVHR = -1;
                else
                    CO_EchoSVHR = ((V_LVdiast_Data - ...    % Calculate SV * HR
                        V_LVsyst_Data) * ...                %  est of cardiac output
                        HR_Echo_Data) / 1000;               %  for 2D echo vols
                    CO_EchoTeich = -1;
                end
                CO_RHCEcho_Max = ...                        % Max of Thermo
                    max([CO_Thermo_Data CO_Fick_Data ...    %  Fick, Teich, SV*HR
                    CO_EchoTeich CO_EchoSVHR ...            %  and EchoD as Norm CO
                    CO_EchoD_Data]);
                % First RHC cardiac output residual
                if (CO_Thermo_Data ~=-1)
                    CO_RHCRes = abs(CO_RHCSim - ...         % Thermo is first choice
                        CO_Thermo_Data) / CO_RHCEcho_Max;
                else
                    CO_RHCRes = abs(CO_RHCSim - ...         % If no thermo then Fick
                        CO_Fick_Data) / CO_RHCEcho_Max;
                end
                % Next Echo cardiac output residual
                if (V_LVsyst_Data == -1 && ...
                    CO_EchoD_Data == -1)
                    CO_EchoRes = abs(CO_EchoSim - ...
                        CO_EchoTeich) / CO_RHCEcho_Max;
                    V_LVsyst_Res = ...
                        abs(V_LVsyst_EchoSim - ...
                        V_LVsyst_Teich) / V_LVdiast_Teich;
                    V_LVdiast_Res = ...
                        abs(V_LVdiast_EchoSim - ...
                        V_LVdiast_Teich) / V_LVdiast_Teich;
                elseif (V_LVsyst_Data == -1 && ...
                        CO_EchoD_Data ~= -1)
                    CO_EchoDRes = abs(CO_EchoSim - ...
                        CO_EchoSVHR) / CO_RHCEcho_Max;
                    CO_EchoTeichRes = abs(CO_EchoSim - ...
                        CO_EchoTeich) / CO_RHCEcho_Max;
                    CO_EchoRes = (CO_EchoTeichRes + ...
                        CO_EchoDRes) / 2;
                    V_LVsyst_Res = ...
                        abs(V_LVsyst_EchoSim - ...
                        V_LVsyst_Teich) / V_LVdiast_Teich;
                    V_LVdiast_Res = ...
                        abs(V_LVdiast_EchoSim - ...
                        V_LVdiast_Teich) / V_LVdiast_Teich;
                else
                    CO_EchoSVHRRes = abs(CO_EchoSim - ...
                        CO_EchoSVHR) / CO_RHCEcho_Max;
                    CO_EchoDRes = abs(CO_EchoSim - ...
                        CO_EchoD_Data) / CO_RHCEcho_Max;
                    CO_EchoRes = (CO_EchoSVHRRes + ...
                        CO_EchoDRes) / 2;
                    V_LVsyst_Res = ...
                        abs(V_LVsyst_EchoSim - ...
                        V_LVsyst_Data) / V_LVdiast_Data;
                    V_LVdiast_Res = ...
                        abs(V_LVdiast_EchoSim - ...
                        V_LVdiast_Data) / V_LVdiast_Data;
                end
                % Now sum up all the residuals and average
                Res = (P_RVsyst_Res + P_RVdiast_Res + ...
                    P_PAsyst_Res + P_PAdiast_Res + ...
                    P_AOsyst_Res + P_AOdiast_Res + ...
                    P_PCWave_Res + CO_RHCRes + ...
                    CO_EchoRes + V_LVsyst_Res + ...
                    V_LVdiast_Res) / 11;
            else
                
                Res = 10;
                
            end
            
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
            % Calculating timespans to reach steady state
            TSpan_SS = [0 NumBeats_SS * period];
            % Solve over the steady state time span with ode15s
            [T_RHC_Out,X_RHC_Out] = ...
                ode15s(@dXdT_SmithRed4,TSpan_SS,X0, ...
                [],DriverP_Struct,CVParam_Struct);
           
            %  Now run Smith et al. model to get intermediate pressures
            Num_TOut_RHC = size(T_RHC_Out,1);    % Number of time points
            P_RV_RHC = zeros(Num_TOut_RHC,1); % Preallocate matrices
            P_AO_RHC = zeros(Num_TOut_RHC,1);
            P_PA_RHC = zeros(Num_TOut_RHC,1);
            P_PU_RHC = zeros(Num_TOut_RHC,1);
            % Use output state variable vector to get pressures
            for i = 1:Num_TOut_RHC
                VarOut = dXdT_Smith(T_RHC_Out(i), ...
                    X_RHC_Out(i,:),DriverP_Struct,CVParam_Struct,1);
                P_RV_RHC(i) = VarOut(2);
                P_AO_RHC(i) = VarOut(3);
                P_PA_RHC(i) = VarOut(5);
                P_PU_RHC(i) = VarOut(6);
            end
            
            % NOW CHECK THAT THE ODES WERE INTEGRATED OVER THE FULL TIME SPAN
            T_SSRHC = NumBeats_SS * (60/HR_RHC_Data);
            if (T_RHC_Out(end) == T_SSRHC)
            
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
                    (NumBeats_Sim - NumBeats_ResPlot) * period;
                tIndStart_RHCResCalc = ...
                    find((T_RHC_Out >= TimeStart_RHCResCalc),1,'first');
                % Looking for systolic and diastolic values in simulation results
                for i = tIndStart_RHCResCalc:Num_TOut_RHCSim
                    P_RVsyst_RHCSim = max(P_RVsyst_RHCSim,P_RV_RHC(i));
                    P_RVdiast_RHCSim = min(P_RVdiast_RHCSim,P_RV_RHC(i));
                    P_PAsyst_RHCSim = max(P_PAsyst_RHCSim,P_PA_RHC(i));
                    P_PAdiast_RHCSim = min(P_PAdiast_RHCSim,P_PA_RHC(i));
                    P_AOsyst_RHCSim = max(P_AOsyst_RHCSim,P_AO_RHC(i));
                    P_AOdiast_RHCSim = min(P_AOdiast_RHCSim,P_AO_RHC(i));
                    P_PCWave_RHCSum = P_PCWave_RHCSum + P_PU_RHC(i);
                    V_LVsyst_RHCSim = min(V_LVsyst_RHCSim,X_RHC_Out(i,1));
                    V_LVdiast_RHCSim = max(V_LVdiast_RHCSim,X_RHC_Out(i,1));
                end
                P_PCWave_RHCSim = P_PCWave_RHCSum / ...
                    (Num_TOut_RHC - tIndStart_RHCResCalc + 1);
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
            else
                
                Res = 10;
                
            end
            
        catch
            
            Res = 10;
            
        end
        
    end
    
end
