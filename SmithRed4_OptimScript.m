% ***********************************************************************************
%         R E D U C E D   V E R S I O N   O F   S M I T H   
%                 C A R D I O V A S C U L A R   S Y S T E M S   M O D E L
% ***********************************************************************************
%
%   This optimization script takes right heart catheter data (RHC) with or without
%   echocardiogram (Echo) data and then simulates a reduced cardiovascular systems 
%   model based on the Smith el al. model (Med Eng Phys 26:131, 2004) that has been
%   matched to this data. The data is drawn from the Cardiovascular Health 
%   Improvement Project (CHIP) data repository at the University of Michigan but any 
%   retrospective or prosepctive dataset with RHC and Echo data can be used. The set 
%   of equations in the original Smith model is a set of differential - algebraic 
%   equations (DAEs) that accounts for ventricular-ventricular interaction (VVI) 
%   however here because of the sparse nature of the clinical data the model has 
%   been reduced to allow the reduced number of parameters to be identified. This
%   results in a reduced form of the Smith model where the VVI, valve inertances, 
%   all dead space/zero pressure volumes, pericardium and thoracic chambers have 
%   been removed.
%
%   Model originally created on     17  January 2016
%   Model last modfied on           19 November 2020
%
%   Developed by        Brian Carlson
%                       Physiological Systems Dynamics Laboratory
%                       Department of Molecular and Integrative Physiology
%                       University of Michigan
%
% ***********************************************************************************
%  Start of                 R E D U C E D   S M I T H   M O D E L   S I M / O P T
% ***********************************************************************************

%% **********************************************************************************
%  Flag Data for            R E D U C E D   S M I T H   M O D E L   S I M / O P T
% ***********************************************************************************

    tic
    warning('off','all')
    
    % Change this to run different patients in the study
    RunPtNum = 27;                                           % Patient number to run
    
    Delim_In = '\t';                                         % Delimiter in file
    Num_HeadLns = 30;                                        % Number of header rows
    FileName2Load = ...                                      % Construct input data
    ['MP' num2str(RunPtNum) '_InputData.txt'];               %  filename
    DataStruct = ...                                         % Load data and header
        importdata(FileName2Load,Delim_In,Num_HeadLns);      %  into structure
    
    if (size(DataStruct.data,1) == 9)                        % Inserting row of NaN
        NaN_Pad = [];                                        %  when actual data
        NumCols = size(DataStruct.data,2);                   %  in file is only
        for i = 1:NumCols                                    %  nine rows. We need
            NaN_Pad = [NaN_Pad NaN];                         %  this because a tenth
        end                                                  %  row is read below
        DataStruct.data = [DataStruct.data; NaN_Pad];
    end
    
    FlagData_Vect = DataStruct.data(1:3,1);                 % Flags controlling code
    PatData_Vect = DataStruct.data(1:4,2);                  % General patient data
    RHCData_Vect = DataStruct.data(1:10,3);                 % RHC data
    if (FlagData_Vect(1) == 1)                              % RHC+Echo data present
        EchoData_Vect = DataStruct.data(1:6,4);             % Echo data
        SimSpecs_Vect = DataStruct.data(1:2,5);             % Hand tune sim beats
    else                                                    % RHC data only present
        SimSpecs_Vect = DataStruct.data(1:2,4);             % Hand tune sim beats
    end
    
    % First we determine how to run this code: check for TTE and RHC data or RHC 
    %  only, run with previously optimized parameter values or parameter values
    %  loaded below, run simulation purely based on nominal parameters
    RHCEcho_Flag = FlagData_Vect(1);            % RHC + Echo = 1, RHC only = 2
    OptimBest_Flag = FlagData_Vect(2);          % 1 = Optim p Run 0 = adjustable run
    NomParam_Flag = FlagData_Vect(3);           % Sim/plot nmnl params = 1, not = 0
    % Save all flag data in a structure to be passed to functions
    FlagData_Values = {RHCEcho_Flag OptimBest_Flag};
    FlagData_Fields = {'RHCEcho_Flag' 'OptimBest_Flag'};
    FlagData_Struct = cell2struct(FlagData_Values, ...
                FlagData_Fields,2);
   
            
%% **********************************************************************************
%  Patient Data for         R E D U C E D   S M I T H   M O D E L   S I M / O P T
% ***********************************************************************************            
           
    % Get general data of study patient number, weight, height and sex
    DID_Num = ...                               % Deidentified subject number
        ['MP' num2str(PatData_Vect(1))];                     
    BW = PatData_Vect(2);                       % Body weight (kg)
    Hgt = PatData_Vect(3);                      % Height (cm)
    if (PatData_Vect(4) == 1)                   % Sex (M or F)
        Sex = 'M';  
    else
        Sex = 'F';
    end
    % Save all patient data into a structure to be passed to functions
    PatData_Values = {DID_Num BW Hgt Sex};
    PatData_Fields = {'DID_Num' 'BW' 'Hgt' 'Sex'};
    PatData_Struct = cell2struct(PatData_Values, ...
        PatData_Fields,2);

    % Now get the RHC data which will be present for all patients
    P_RVsyst = RHCData_Vect(1);                 % Systolic RV pressure (mmHg) 
    P_RVdiast = RHCData_Vect(2);                % Diastolic RV pressure (mmHg)
    P_PAsyst = RHCData_Vect(3);                 % Systolic pulm art press (mmHg)
    P_PAdiast = RHCData_Vect(4);                % Diast pulm art press (mmHg)
    P_PCWave = RHCData_Vect(5);                 % Average pulm wedge press (mmHg)
    P_AOsyst = RHCData_Vect(6);                 % Systolic aortic press (mmHg)
    P_AOdiast = RHCData_Vect(7);                % Diastolic aortic press (mmHg)
    HR_RHC = RHCData_Vect(8);                   % Average heart rate (beats/min)
    CO_Fick = RHCData_Vect(9);                  % Cardiac output Fick (L/min)
    if (isnan(RHCData_Vect(10)))
        CO_Thermo = -1;                         % No CO thermo in RHC record
    else
        CO_Thermo = RHCData_Vect(10);           % Cardiac outpt thermodil (L/min)
    end
    % Save all RHC data into a structure to be passed to functions
    RHCData_Values = {P_RVsyst P_RVdiast P_PAsyst P_PAdiast ....
        P_PCWave P_AOsyst P_AOdiast HR_RHC CO_Fick CO_Thermo};
    RHCData_Fields = {'P_RVsyst' 'P_RVdiast' 'P_PAsyst' 'P_PAdiast' ....
        'P_PCWave' 'P_AOsyst' 'P_AOdiast' 'HR_RHC' 'CO_Fick' 'CO_Thermo'};
    RHCData_Struct = cell2struct(RHCData_Values, ...
        RHCData_Fields,2);

    % Get the echo data if present in the input file
    if (RHCEcho_Flag == 1)
        ID_LVsyst = EchoData_Vect(1);           % Systolic LV inner diam (mm)
        ID_LVdiast = EchoData_Vect(2);          % Diastolic LV inner diam (mm)
        HR_Echo	= EchoData_Vect(3);             % Average heart rate (beats/min)
        if (isnan(EchoData_Vect(4)))
            CO_EchoD = -1;                      % No Echo-Dop cardiac output
        else
            CO_EchoD = EchoData_Vect(4);        % Cardiac output Echo-Dop (L/min)
        end
        if (isnan(EchoData_Vect(5)))            
            V_LVsyst = -1;                      % 2D echo volumes
            V_LVdiast = -1;                     %  were not calculated
        else
            V_LVsyst = EchoData_Vect(5);        % Systolic LV volume (mL)
            V_LVdiast = EchoData_Vect(6);       % Diastolic LV volume (mL)
        end
        % Save all Echo data into a structure to be passed to functions
        EchoData_Values = {ID_LVsyst ID_LVdiast HR_Echo ...
            CO_EchoD V_LVsyst V_LVdiast};
        EchoData_Fields = {'ID_LVsyst' 'ID_LVdiast' 'HR_Echo' ...
            'CO_EchoD' 'V_LVsyst' 'V_LVdiast'};
        EchoData_Struct = cell2struct(EchoData_Values, ...
            EchoData_Fields,2);
    end
        
    

%% **********************************************************************************
%  Fixed Params of          R E D U C E D   S M I T H   M O D E L   S I M / O P T
% ***********************************************************************************
    
    % NORMAL MODEL PARAMETERS (* MEANS TAKEN DIRECTLY FROM SMITH ET AL.)
    % This set of parameters give roughly 120/80 mmHg systemic pressure, 
    %  20/9 mmHg pulmonary artery pressure, left ventricular diastolic 
    %  volume of 85 mL, left ventricular stroke volume of 57mL (for an
    %  ejection fraction of 67%) and cardiac output of 4.6 L/min. Parameters
    %  that have a multiplicative factor indicate the relative change from
    %  the full Smith et al. model.
    if (OptimBest_Flag == 0)
            
        % Elastance function driver parameters
        A = 1;                                      % Elastance funct param (uls)   * 
        % Left ventricle free wall parameters
        E_es_lvf = 2.8798 * 1.50;                   % LV free wall elast (mmHg/mL)  *
        P_0_lvf = 0.1203;                           % LV ED pressure param (mmHg)   *
        lambda_lvf = 0.033 * 0.70; % * 0.325;       % LV ED pressure param (1/mL)
        % Right ventricle free wall parameters
        E_es_rvf = 0.585 * 1.20; % * 0.85;          % RV free wall elstnce (mmHg/mL)
        P_0_rvf = 0.2157;                           % RV ED pressure param (mmHg)   *
        lambda_rvf = 0.023 * 0.70; % * 0.325;       % RV ED pressure param (1/mL)
        % Pulmonary artery and vein parameters
        E_es_pa = 0.369 * 0.70;                     % Pulm arter elstnce (mmHg/mL)  *
        E_es_pu = 0.0073;                           % Pulm ven elastance (mmHg/mL)  *
        R_pul = 0.1552 * 0.85; % * 0.65;            % Pulm vasc resist (mmHg*s/mL)
        % Aortic and vena cava parameters
        E_es_sa = 0.6913 * 1.30; % * 1.20;          % Syst arter elstnce (mmHg/mL)
        E_es_sv = 0.0059;                           % Syst venous elstnce (mmHg/mL) *
        R_sys = 1.0889 * 1.18; % * 1.11;            % Syst vasc resist (mmHg*s/mL)
        % Heart valve paramenters
        R_mt = 0.0158;                              % Mitral vlv resist (mmHg*s/mL) *
        R_av = 0.018;                               % Aortic vlv resist (mmHg*s/mL) *
        R_tc = 0.0237;                              % Tricspd vlv resist (mmHg*s/mL)*
        R_pv = 0.0055;                              % Pulmon vlv resist (mmHg*s/mL) *
        % Heart failure param
        SVFact = 1.00;                              % Stress blood vol factor (uls)
        
        % Save all parameters in a structure to pass
        CVParam_Values = {A E_es_lvf P_0_lvf lambda_lvf E_es_rvf ...
            P_0_rvf lambda_rvf E_es_pa E_es_pu R_pul E_es_sa ...
            E_es_sv R_sys R_mt R_av R_tc R_pv SVFact};

        CVParam_Fields = {'A' 'E_es_lvf' 'P_0_lvf' 'lambda_lvf' ...
            'E_es_rvf' 'P_0_rvf' 'lambda_rvf' 'E_es_pa' ...
            'E_es_pu' 'R_pul' 'E_es_sa' 'E_es_sv' 'R_sys' ...
            'R_mt' 'R_av' 'R_tc' 'R_pv' 'SVFact'};
        CVParam_Struct = cell2struct(CVParam_Values,CVParam_Fields,2);

        % In this case we are not simulating the best optimization or
        %  doing a hand tune of model to data but are running a simulation
        %  using the nominal parameters
        if (NomParam_Flag == 1)

            % Set these to normal cardiovascular function values (see above)
            %  since they are not recalculated in the nominal parameter
            %  calculation function but just set to normal values
            P_0_lvf = 0.1203;                   % LV ED pressure param (mmHg)
            P_0_rvf = 0.1203;                   % RV ED pressure param (mmHg)
            SVFact = 1.0000;                    % Stress blood vol factor (uls)
            % Now load them into the structure so they can be 
            %  used in the nominal parameter calculation function
            CVParam_Struct.P_0_lvf = P_0_lvf;
            CVParam_Struct.P_0_lvf = P_0_lvf;
            CVParam_Struct.SVFact = SVFact;

            % OPTIMIZATION IS BEING RUN NEEDING NOMINAL PARAMETER CALCULATIONS
            NomParam_Struct = NomParam_Calc(PatData_Struct, ...
                CVParam_Struct,RHCData_Struct,EchoData_Struct);
            
            % Overwrite the CVParam structure with nominal parameter values
            %  for simulation and visulaization
            CVParam_Struct.E_es_lvf = NomParam_Struct.Eeslvf_Nom;
            CVParam_Struct.P_0_lvf = NomParam_Struct.P0lvf_Nom;
            CVParam_Struct.lambda_lvf = NomParam_Struct.lambdalvf_Nom;
            CVParam_Struct.E_es_rvf = NomParam_Struct.Eesrvf_Nom;
            CVParam_Struct.P_0_rvf = NomParam_Struct.P0rvf_Nom;
            CVParam_Struct.lambda_rvf = NomParam_Struct.lambdarvf_Nom;
            CVParam_Struct.E_es_pa = NomParam_Struct.Eespa_Nom;
            CVParam_Struct.E_es_pu = NomParam_Struct.Eespu_Nom;
            CVParam_Struct.R_pul = NomParam_Struct.Rpul_Nom;
            CVParam_Struct.E_es_sa = NomParam_Struct.Eessa_Nom;
            CVParam_Struct.E_es_sv = NomParam_Struct.Eessv_Nom;
            CVParam_Struct.R_sys = NomParam_Struct.Rsys_Nom;
            CVParam_Struct.R_mt = NomParam_Struct.Rmt_Nom;
            CVParam_Struct.R_av = NomParam_Struct.Rav_Nom;
            CVParam_Struct.R_tc = NomParam_Struct.Rtc_Nom;
            CVParam_Struct.R_pv = NomParam_Struct.Rpv_Nom;

        end
                
        
    else
        
        % LOAD THE OPTIM PARAMS FROM AN EARLIER RUN
        Delim_In = '\t';
        Num_HeadLns = 22;
        BestParams2Load = ['MP' num2str(RunPtNum) '_OptimBestParams.txt'];
        BestParams_Struct = importdata(BestParams2Load,Delim_In,Num_HeadLns);
        % Elastance function driver parameters
        A = BestParams_Struct.data(1);          % Elastance function param (uls)
        % Left ventricle free wall parameters
        E_es_lvf = BestParams_Struct.data(2);   % LV free wall elast (mmHg/mL) 
        P_0_lvf = BestParams_Struct.data(3);    % LV ED pressure param (mmHg)
        lambda_lvf = BestParams_Struct.data(4); % LV ED pressure param (1/mL)
        % Right ventricle free wall parameters
        E_es_rvf = BestParams_Struct.data(5);   % RV free wall elast (mmHg/mL) 
        P_0_rvf = BestParams_Struct.data(6);    % RV ED pressure param (mmHg)
        lambda_rvf = BestParams_Struct.data(7); % RV ED pressure param (1/mL)
        % Pulmonary artery and vein parameters
        E_es_pa = BestParams_Struct.data(8);    % Pulm arterial elstnce (mmHg/mL)
        E_es_pu = BestParams_Struct.data(9);    % Pulm venous elastance (mmHg/mL)
        R_pul = BestParams_Struct.data(10);     % Pulm vasc resist (mmHg*s/mL)
        % Aortic and vena cava parameters
        E_es_sa = BestParams_Struct.data(11);   % Syst arterial elstnce (mmHg/mL)
        E_es_sv = BestParams_Struct.data(12);   % Syst venous elastance (mmHg/mL)
        R_sys = BestParams_Struct.data(13);     % Syst vasc resistnce (mmHg*s/mL)
        % Heart valve paramenters
        R_mt = BestParams_Struct.data(14);      % Mitral valve resist (mmHg*s/mL)
        R_av = BestParams_Struct.data(15);      % Aortic valve resist (mmHg*s/mL)
        R_tc = BestParams_Struct.data(16);      % Tricspd vlv resist (mmHg*s/mL)
        R_pv = BestParams_Struct.data(17);      % Pulmon vlv resist (mmHg*s/mL)
        % Heart failure param
        SVFact = BestParams_Struct.data(18);    % Stress blood vol factor (uls)

        % Save all parameters in a structure to pass
        CVParam_Values = {A E_es_lvf P_0_lvf lambda_lvf ...
            E_es_rvf P_0_rvf lambda_rvf E_es_pa E_es_pu R_pul ...
            E_es_sa E_es_sv R_sys R_mt R_av R_tc R_pv SVFact};

        CVParam_Fields = {'A' 'E_es_lvf' 'P_0_lvf' 'lambda_lvf' ...
            'E_es_rvf' 'P_0_rvf' 'lambda_rvf' 'E_es_pa' ...
            'E_es_pu' 'R_pul' 'E_es_sa' 'E_es_sv' 'R_sys' ...
            'R_mt' 'R_av' 'R_tc' 'R_pv' 'SVFact'};
        CVParam_Struct = cell2struct(CVParam_Values,CVParam_Fields,2);

    end
            

%% **********************************************************************************
%  Opt/Sim Params of        R E D U C E D   S M I T H   M O D E L   S I M / O P T
% ***********************************************************************************

    % Simulation parameters
    NumBeats_SS = SimSpecs_Vect(1);                      % Num beats 2 stdy state
    NumBeats_ResPlot = SimSpecs_Vect(2);                 % Num beats to plot

    % Building the simulation parameter structure
    SimOptParam_Values = {NumBeats_SS NumBeats_ResPlot};
    SimOptParam_Fields = {'NumBeats_SS' 'NumBeats_ResPlot'};
    SimOptParam_Struct = cell2struct(SimOptParam_Values, ...
        SimOptParam_Fields,2);
    
    % Now put all structures into a single structure to pass   
    if (RHCEcho_Flag == 1)
        AllStruct_Values = {FlagData_Struct PatData_Struct ...
            RHCData_Struct EchoData_Struct ...
            CVParam_Struct SimOptParam_Struct};
        AllStruct_Fields = {'FlagData_Struct' 'PatData_Struct' ...
            'RHCData_Struct' 'EchoData_Struct' ...
            'CVParam_Struct' 'SimOptParam_Struct'};
        AllStruct_Struct = cell2struct(AllStruct_Values, ...
        AllStruct_Fields,2);
    else
        AllStruct_Values = {FlagData_Struct PatData_Struct ...
            RHCData_Struct CVParam_Struct SimOptParam_Struct};
        AllStruct_Fields = {'FlagData_Struct' 'PatData_Struct' ...
            'RHCData_Struct' 'CVParam_Struct' 'SimOptParam_Struct'};
        AllStruct_Struct = cell2struct(AllStruct_Values, ...
        AllStruct_Fields,2);
    end
    
    
%% **********************************************************************************
%  Run Simulation of        R E D U C E D   S M I T H   M O D E L   S I M / O P T
% ***********************************************************************************
        
    % Calculate total blood volume based on height, weight and sex
    %  which is the same for both RHC and Echo simulations
    %  This expression is from Nadler et al. Surgery 51:224,1962.
    if (Sex == 'M')
        TotBV = ((0.3669 * (Hgt/100)^3) + ...
            (0.03219 * BW) + 0.6041) * 1000;
    else
        TotBV = ((0.3561 * (Hgt/100)^3) + ...
            (0.03308 * BW) + 0.1833) * 1000;
    end
    % The original Smith model accounted for the total stressed blood volume as
    %  a fraction of the total blood volume. Assuming they were simulating a typical 
    %  5000 mL total blood volume they included only 1500 mL (or 30%) as the 
    %  stressed volume therefore we will multiply our calculated TotBV value by 0.3 
    %  to yield the stressed blood volume which the equations keep track of in each
    %  compartment at each time step. To account for extra recruited stressed volume 
    %  that may be present in heart disease the 30% stressed blood volume can be 
    %  altered by changing SVFact. The cardiovascular dynamics are greatly affected 
    %  by stressed volume fraction so in this study we held the total stressed 
    %  volume constant at 30% of total blood volume by fixing SVFact at 1.00.
    SVFact = CVParam_Struct.SVFact;
    CircBV = SVFact * 0.30 * TotBV;
    % Setting state variable initial conditions 
    %  Note that initial volume division is scaled as a fraction
    %  of stressed blood volume calculated earlier and the 
    %  initial guess of flow is total blood volume in one minute
    V_lv0 = (94.6812/1500) * CircBV;
    V_rv0 = (90.7302/1500) * CircBV;
    V_pa0 = (43.0123/1500) * CircBV;
    V_pu0 = (808.458/1500) * CircBV;
    V_sa0 = (133.338/1500) * CircBV;
    V_sv0 = (329.780/1500) * CircBV;
    % Put into vector to pass to ode15s
    X0(1) = V_lv0;
    X0(2) = V_rv0;
    X0(3) = V_pa0;
    X0(4) = V_pu0;
    X0(5) = V_sa0;
    X0(6) = V_sv0;

    if (RHCEcho_Flag == 1)                      % With RHC and Echo data

        % RUN RHC SIMULATION FIRST
        % Build driver function parameter structure
        HR = HR_RHC;                            % RHC heart rate (beats/min)
        period_RHC = 60/HR_RHC;                 % Period of heart beat (s)
        B_RHC = HR_RHC;                         % Elastance fctn param (1/s^2)
        C_RHC = period_RHC/2;                   % Elastance fctn param (s)
        DriverP_Values = {HR period_RHC B_RHC C_RHC};
        DriverP_Fields = {'HR' 'period' 'B' 'C'};
        DriverP_Struct = cell2struct(DriverP_Values, ...
            DriverP_Fields,2);
        % Calculating timespans to reach steady state and then for simulation
        TSpan_SSRHC = [0 NumBeats_SS * period_RHC];

        % Solve over the time span with ode15s
        [T_Out_RHC,X_Out_RHC] = ode15s(@dXdT_SmithRed4, ...
            TSpan_SSRHC,X0,[],DriverP_Struct,CVParam_Struct);

        % CAPTURING INTERMEDIATE PRESSURES TO PLOT
        Num_TOut_RHC = size(T_Out_RHC,1); % Number of time points
        P_LV_RHC = zeros(Num_TOut_RHC,1); % Preallocating matrices
        P_RV_RHC = zeros(Num_TOut_RHC,1);
        P_AO_RHC = zeros(Num_TOut_RHC,1);
        P_VC_RHC = zeros(Num_TOut_RHC,1);
        P_PA_RHC = zeros(Num_TOut_RHC,1);
        P_PU_RHC = zeros(Num_TOut_RHC,1);
        % RERUNNING MODEL TO GET INTERMEDIATE PRESSURES
        for i = 1:Num_TOut_RHC
            VarOut = dXdT_SmithRed4(T_Out_RHC(i), ...
                X_Out_RHC(i,:),DriverP_Struct, ...
                CVParam_Struct,1);
            P_LV_RHC(i) = VarOut(1);
            P_RV_RHC(i) = VarOut(2);
            P_AO_RHC(i) = VarOut(3);
            P_VC_RHC(i) = VarOut(4);
            P_PA_RHC(i) = VarOut(5);
            P_PU_RHC(i) = VarOut(6);
        end

        % RUN ECHO SIMULATION NEXT
        % Build driver function parameter structure
        HR = HR_Echo;                           % RHC heart rate (beats/min)
        period_Echo = 60/HR_Echo;               % Period of heart beat (s)
        B_Echo = HR_Echo;                       % Elastance fctn param (1/s^2)
        C_Echo = period_Echo/2;                 % Elastance fctn param (s)
        DriverP_Values = {HR period_Echo B_Echo C_Echo};
        DriverP_Fields = {'HR' 'period' 'B' 'C'};
        DriverP_Struct = cell2struct(DriverP_Values, ...
            DriverP_Fields,2);
        % Calculating timespans to reach steady state and then for simulation
        TSpan_SSEcho = [0 NumBeats_SS * period_Echo];

        % Solve over the time span with ode15s
        [T_Out_Echo,X_Out_Echo] = ode15s(@dXdT_SmithRed4, ...
            TSpan_SSEcho,X0,[],DriverP_Struct,CVParam_Struct);

        % CAPTURING THE LEFT AND RIGHT VENTRICULAR PRESSURES
        Num_TOut_Echo = ...                     % Number of time points
            size(T_Out_Echo,1);
        P_LV_Echo = zeros(Num_TOut_Echo,1);     % Preallocating matrices
        P_RV_Echo = zeros(Num_TOut_Echo,1);
        % RERUNNING THE SIMULATION TO EXTRACT THE PRESSURES
        for i = 1:Num_TOut_Echo
            VarOut = dXdT_SmithRed4(T_Out_Echo(i), ...
                X_Out_Echo(i,:),DriverP_Struct, ...
                CVParam_Struct,1);
            P_LV_Echo(i) = VarOut(1);
            P_RV_Echo(i) = VarOut(2);
        end

    else                                        % RHC data only

        % RUN RHC SIMULATION ONLY
        % Build driver function parameter structure
        HR = HR_RHC;                            % RHC heart rate (beats/min)
        period_RHC = 60/HR_RHC;                 % Period of heart beat (s)
        B_RHC = HR_RHC;                         % Elastance fctn param (1/s^2)
        C_RHC = period_RHC/2;                   % Elastance fctn param (s)
        DriverP_Values = {HR period_RHC B_RHC C_RHC};
        DriverP_Fields = {'HR' 'period' 'B' 'C'};
        DriverP_Struct = cell2struct(DriverP_Values, ...
            DriverP_Fields,2);
        %Calculating timespans to reach steady state and then for simulation
        TSpan_SS = [0 NumBeats_SS * period_RHC];

        % Solve over the time span with ode15s
        [T_Out_RHC,X_Out_RHC] = ode15s(@dXdT_SmithRed4, ...
            TSpan_SS,X0,[],DriverP_Struct,CVParam_Struct);

        % CAPTURING INTERMEDIATE PRESSURES TO PLOT
        Num_TOut_RHC = size(T_Out_RHC,1); % Number of time points
        P_LV_RHC = zeros(Num_TOut_RHC,1); % Preallocating matrices
        P_RV_RHC = zeros(Num_TOut_RHC,1);
        P_AO_RHC = zeros(Num_TOut_RHC,1);
        P_VC_RHC = zeros(Num_TOut_RHC,1);
        P_PA_RHC = zeros(Num_TOut_RHC,1);
        P_PU_RHC = zeros(Num_TOut_RHC,1);
        % RERUNNING MODEL TO GET INTERMEDIATE PRESSURES
        for i = 1:Num_TOut_RHC
            VarOut = dXdT_SmithRed4(T_Out_RHC(i), ...
                X_Out_RHC(i,:),DriverP_Struct, ...
                CVParam_Struct,1);
            P_LV_RHC(i) = VarOut(1);
            P_RV_RHC(i) = VarOut(2);
            P_AO_RHC(i) = VarOut(3);
            P_VC_RHC(i) = VarOut(4);
            P_PA_RHC(i) = VarOut(5);
            P_PU_RHC(i) = VarOut(6);
        end

    end
        

%% **********************************************************************************
%  Plot Sim or Opt of   S M I T H   C A R D I O V A S C   M O D E L   S I M / O P T
% ***********************************************************************************
        
    if (RHCEcho_Flag == 1)

        % PATIENT SPECIFIC MODEL RESULTS - ONLY RHC + TTE
        % Finding the index and the time associated with the last
        %  few beats to plot which were used for the residual calculation
        NumBeat_Start = NumBeats_SS - NumBeats_ResPlot; % Starting beat
        tStartRHC_Ind = find(T_Out_RHC >= ...           
            (NumBeat_Start * period_RHC),1,'first');
        tStartRHC = T_Out_RHC(tStartRHC_Ind);
        tStartEcho_Ind = ...
            find(T_Out_Echo >= (NumBeat_Start * period_Echo),1,'first');
        tStartEcho = T_Out_Echo(tStartEcho_Ind);
        T_OutRHC = T_Out_RHC(tStartRHC_Ind:end) - tStartRHC;
        T_OutEcho = T_Out_Echo(tStartEcho_Ind:end) - tStartEcho;
        P_RVRHCSim = P_RV_RHC(tStartRHC_Ind:end);
        P_RVEchoSim = P_RV_Echo(tStartEcho_Ind:end);
        P_AOSim = P_AO_RHC(tStartRHC_Ind:end);
        P_PASim = P_PA_RHC(tStartRHC_Ind:end);
        P_PUSim = P_PU_RHC(tStartRHC_Ind:end);
        V_LVRHCSim = X_Out_RHC(tStartRHC_Ind:end,1);
        V_RVRHCSim = X_Out_RHC(tStartRHC_Ind:end,2);
        V_LVED_RHC = max(V_LVRHCSim);
        V_LVES_RHC = min(V_LVRHCSim);    
        CO_RHCSim = ((V_LVED_RHC - V_LVES_RHC) * HR_RHC) / 1000;
        V_LVEchoSim = X_Out_Echo(tStartEcho_Ind:end,1);
        V_RVEchoSim = X_Out_Echo(tStartEcho_Ind:end,2);
        V_LVED_Echo = max(V_LVEchoSim);
        V_LVES_Echo = min(V_LVEchoSim);
        CO_EchoSim = ((V_LVED_Echo - V_LVES_Echo) * HR_Echo) / 1000;
        P_LVRHCSim = P_LV_RHC(tStartRHC_Ind:end);
        P_LVEchoSim = P_LV_Echo(tStartEcho_Ind:end);
        P_VCSim = P_VC_RHC(tStartRHC_Ind:end);

        SimPlot_Values = {T_OutRHC T_OutEcho P_RVRHCSim ...
            P_RVEchoSim P_AOSim P_PASim P_PUSim CO_RHCSim ...
            CO_EchoSim V_LVEchoSim V_RVEchoSim P_LVRHCSim ...
            P_LVEchoSim P_VCSim};
        SimPlot_Fields = {'T_OutRHC' 'T_OutEcho' 'P_RVRHCSim' ...
            'P_RVEchoSim' 'P_AOSim' 'P_PASim' 'P_PUSim' ...
            'CO_RHCSim' 'CO_EchoSim' 'V_LVEchoSim' 'V_RVEchoSim' ...
            'P_LVRHCSim' 'P_LVEchoSim' 'P_VCSim'};
        SimPlot_Struct = cell2struct(SimPlot_Values, ...
            SimPlot_Fields,2);

        SixPanel_Figure(AllStruct_Struct,SimPlot_Struct)

    else

        NumBeat_Start = NumBeats_SS - NumBeats_ResPlot;
        tStartRHC_Ind = ...
            find(T_Out_RHC >= (NumBeat_Start * period_RHC),1,'first');
        T_Out = T_Out_RHC(tStartRHC_Ind:end);
        P_RVSim = P_RV_RHC(tStartRHC_Ind:end);
        P_AOSim = P_AO_RHC(tStartRHC_Ind:end);
        P_PASim = P_PA_RHC(tStartRHC_Ind:end);
        P_PUSim = P_PU_RHC(tStartRHC_Ind:end);
        V_LVSim = X_Out_RHC(tStartRHC_Ind:end,1);
        V_RVSim = X_Out_RHC(tStartRHC_Ind:end,2);
        P_LVSim = P_LV_RHC(tStartRHC_Ind:end);
        P_VCSim = P_VC_RHC(tStartRHC_Ind:end);

        SimPlot_Values = {T_Out P_RVSim P_AOSim P_PASim ...
            P_PUSim CO_RHCSim V_LVSim V_RVSim P_LVSim P_VCSim};
        SimPlot_Fields = {'T_Out' 'P_RVSim' 'P_AOSim' 'P_PASim' ...
            'P_PUSim' 'CO_RHCSim' 'V_LVSim' 'V_RVSim' ...
            'P_LVSim' 'P_VCSim'};
        SimPlot_Struct = cell2struct(SimPlot_Values, ...
            SimPlot_Fields,2);

        SixPanel_Figure(AllStruct_Struct,SimPlot_Struct)

    end
    

    toc