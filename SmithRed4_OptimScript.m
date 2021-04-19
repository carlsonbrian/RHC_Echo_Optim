% ***********************************************************************************
%         R E D U C E D   V E R S I O N   O F   S M I T H   
%                 C A R D I O V A S C U L A R   S Y S T E M S   M O D E L
% ***********************************************************************************
%
%   This optimization script takes right heart catheter data (RHC) with or without
%   echocardiogram (Echo) data and then optimizes a cardiovascular systems model 
%   based on the Smith el al. model (Med Eng Phys 26:131, 2004) to match this data. 
%   The data is currently being drawn from the Deidentified Clinical Data Repository
%   (DCDR) at the University of Washington or the Cardiovascular Health Improvement 
%   Project (CHIP) data repository at the University of Michigan but any retro-
%   spective or prosepctive dataset with RHC and Echo data can be used. The set of 
%   equations in the Smith model is a set of differential - algebraic 
%   equations (DAEs) and the dXdT called contains a single expression where the
%   right hand side is equal to 0 not the dX/dt of a state variable. This means 
%   when ode15s is called it loads a singular mass matrix M with ones along the 
%   diagonal except for a zero in the position of the implicit expression. In this 
%   optimization script a reduced form of the Smith model is used where the 
%   ventricular-ventricular interaction, valve inertances, all dead space/zero 
%   pressure volumes, pericardium and thoracic chambers have been removed.
%
%   Model originally created on     17 January 2016
%   Model last modfied on            6    July 2020
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
    RunPtNum = 3;                                            % Patient number to run
    
    Delim_In = '\t';                                         % Delimiter in file
    Num_HeadLns = 43;                                        % Number of header rows
    FileName2Load = ...                                      % Construct input data
    ['BCHF' num2str(RunPtNum) 'Red4_InputData.txt'];         %  filename
    DataStruct = ...                                         % Load data and header
        importdata(FileName2Load,Delim_In,Num_HeadLns);      %  into structure
    
    if (size(DataStruct.data,1) == 9)                        % Inserting row of NaN
        NaN_Pad = [];                                        %  when actual data
        NumCols = size(DataStruct.data,2);                   %  in file is only
        for i = 1:NumCols                                    %  nine rows. We need
            NaN_Pad = [NaN_Pad NaN];                         %  this beacuse a tenth
        end                                                  %  row is read below
        DataStruct.data = [DataStruct.data; NaN_Pad];
    end
    
    FlagData_Vect = DataStruct.data(1:9,1);                 % Flags controlling code
    PatData_Vect = DataStruct.data(1:4,2);                  % General patient data
    RHCData_Vect = DataStruct.data(1:10,3);                 % RHC data
    if (FlagData_Vect(1) == 1)                              % RHC+Echo data present
        EchoData_Vect = DataStruct.data(1:6,4);             % Echo data
        if (FlagData_Vect(2) == 1)
            SimOptSpecs_Vect = DataStruct.data(1:2,5);      % Hand tune sim beats
        else                                                %  and plot beats 
            if (FlagData_Vect(3) == 1)
                SimOptSpecs_Vect = DataStruct.data(1:4,5);  % fmc optim sim, optim
            else                                            %  beats and fmc specs
                SimOptSpecs_Vect = DataStruct.data(1:6,5);  % GA optim sim, optim
            end                                             %  beats and GA specs
            Num_DataRows = size(DataStruct.data,1);         % Find num of adj params
            AdjParam_Vect = ...                             % Load adj parameter
                DataStruct.data(1:Num_DataRows,6);          %  numbers
        end
    else                                                    % RHC data only present
        if (FlagData_Vect(2) == 1)
            SimOptSpecs_Vect = DataStruct.data(1:2,4);      % Hand tune sim beats
        else                                                %  and plot beats 
            if (FlagData_Vect(3) == 1)
                SimOptSpecs_Vect = DataStruct.data(1:4,4);  % fmc optim sim, optim
            else                                            %  beats and fmc specs
                SimOptSpecs_Vect = DataStruct.data(1:6,4);  % GA optim sim, optim
            end                                             %  beats and GA specs
            Num_DataRows = size(DataStruct.data,1);         % Find num of adj params
            AdjParam_Vect = ...                             % Load adj parameter
                DataStruct.data(1:Num_DataRows,5);          %  numbers
        end
    end
    
    % First we determine what we want to do with this code: optimize to data,
    %  compare to Smith model simulation, hand fit, run Smith model params and
    %  save output, or plot out a previous optimized fit. Specific detail on
    %  flag values for each option are given above in the script header
    RHCEcho_Flag = FlagData_Vect(1);            % RHC + Echo = 1, RHC only = 2
    HandTune_Flag = FlagData_Vect(2);           % 1 = Hnd tne or sim, 0 = Opt & sim
    Optim_Flag = FlagData_Vect(3);              % 1 = fmincon, 2 = Genetic alg
    Parallel_Flag = FlagData_Vect(4);           % 1 = Parallel, 0 = Serial optim
    OptimBest_Flag = FlagData_Vect(5);          % 1 = Optim p Run 0 = Opt or Other
    NormParam_Flag = FlagData_Vect(6);          % Run normal params = 1, not = 0
    SaveNorm_Flag = FlagData_Vect(7);           % Save normal run = 1, not = 0
    Comp2Norm_Flag = FlagData_Vect(8);          % Compare to normal = 1, not = 0 
    NomParam_Flag = FlagData_Vect(9);           % Sim/plot nmnl params = 1, not = 0
    % Save all flag data in a structure to be passed to functions
    FlagData_Values = {HandTune_Flag OptimBest_Flag Optim_Flag ...
        Parallel_Flag RHCEcho_Flag NormParam_Flag ...
        Comp2Norm_Flag SaveNorm_Flag};
    FlagData_Fields = {'HandTune_Flag' 'OptimBest_Flag' ...
        'Optim_Flag' 'Parallel_Flag' 'RHCEcho_Flag' ...
        'NormParam_Flag' 'Comp2Norm_Flag' 'SaveNorm_Flag'};
    FlagData_Struct = cell2struct(FlagData_Values, ...
                FlagData_Fields,2);
   
%% **********************************************************************************
%  Patient Data for         R E D U C E D   S M I T H   M O D E L   S I M / O P T
% ***********************************************************************************            
            
    if (NormParam_Flag == 0)
           
        % We are running patient data simulation or optimization so get the 
        % general data of study patient number, weight, height and sex
        DID_Num = ...                               % Deidentified subject number
            ['BCHF' num2str(PatData_Vect(1))];                     
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
        
    else
        
        DID_Num = 'Normal';
        HR_Norm = 80;                               % Average heart rate (beats/min)
        BW = 72.3;                                  % Body weight (kg, ~160 lbs)
        Hgt = 178;                                  % Height (cm, ~5'10")
        Sex = 'M';                                  % Assuming male since TotBV ~ 5L
        
        PatData_Values = {DID_Num HR_Norm BW Hgt Sex};
        PatData_Fields = {'DID_Num' 'HR_Norm' 'BW' 'Hgt' 'Sex'};
        PatData_Struct = cell2struct(PatData_Values, ...
            PatData_Fields,2);
        
    end
    

%% **********************************************************************************
%  Fixed Params of          R E D U C E D   S M I T H   M O D E L   S I M / O P T
% ***********************************************************************************
    
    if (NormParam_Flag == 1)
         
        % NORMAL MODEL PARAMETERS (* MEANS TAKEN DIRECTLY FROM SMITH ET AL.)
        % This set of parameters give roughly 120/80 mmHg systemic pressure, 
        %  20/9 mmHg pulmonary artery pressure, left ventricular diastolic 
        %  volume of 85 mL, left ventricular stroke volume of 57mL (for an
        %  ejection fraction of 67%) and cardiac output of 4.6 L/min. Parameters
        %  that have a multiplicative factor indicate the relative change from
        %  the full Smith et al. model.
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
            'E_es_rvf' 'P_0_rvf' 'lambda_rvf' 'E_es_pa' 'E_es_pu' ...
            'R_pul' 'E_es_sa' 'E_es_sv' 'R_sys' 'R_mt' 'R_av' ...
            'R_tc' 'R_pv' 'SVFact'};
        CVParam_Struct = cell2struct(CVParam_Values,CVParam_Fields,2);
        
    else
        
        if (OptimBest_Flag == 0)
            
            % MODEL PARAMETERS TO ADJUST BY HAND OR CALCULATE NOMINALLY
            % Elastance function driver parameters
            A = 1;                                  % Elastance function param (uls)
            % Left ventricle free wall parameters
            E_es_lvf = 2.246866; % 2.8798;          % LV free wall elast (mmHg/mL) 
            P_0_lvf = 0.1203; % 0.1203;             % LV ED pressure param (mmHg)
            lambda_lvf = 0.033; % 0.033;            % LV ED pressure param (1/mL)
            % Right ventricle free wall parameters
            E_es_rvf = 0.859633; % 0.585;           % RV free wall elast (mmHg/mL) 
            P_0_rvf = 0.2157; % 0.2157;             % RV ED pressure param (mmHg)
            lambda_rvf = 0.023; % 0.023;            % RV ED pressure param (1/mL)
            % Pulmonary artery and vein parameters
            E_es_pa = 0.736206; % 0.369;            % Pulm arterial elstnce (mmHg/mL)
            E_es_pu = 0.087768; % 0.0073;           % Pulm venous elastance (mmHg/mL)
            R_pul = 0.1168237; % 0.1552;            % Pulm vasc resist (mmHg*s/mL)
            % Aortic and vena cava parameters
            E_es_sa = 0.5598488; % 0.6913;          % Syst arterial elstnce (mmHg/mL)
            E_es_sv = 0.00971877; % 0.0059;         % Syst venous elastance (mmHg/mL)
            R_sys = 0.836003 * 0.80; % 1.0889;             % Syst vasc resistnce (mmHg*s/mL)
            % Heart valve paramenters
            R_mt = 0.00479139; % 0.0158;            % Mitral valve resist (mmHg*s/mL)
            R_av = 0.02402423; % 0.018;             % Aortic valve resist (mmHg*s/mL)
            R_tc = 0.02624495; % 0.0237;            % Tricspd vlv resist (mmHg*s/mL)
            R_pv = 0.00847914; % 0.0055;            % Pulmon vlv resist (mmHg*s/mL)
            % Heart failure param
            SVFact = 1.00;                          % Stress blood vol factor (uls)
            
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
            %  doing a hand tune of model to data but are running an
            %  optimization where we need nominal parameter values calculated or
            %  simply plotting out the nominal parameter fit
            if (HandTune_Flag == 0)
                
                % Set these to normal cardiovascular function values (see above)
                %  since they are not recalculated in the nominal parameter
                %  calculation function but just set to normal values
                P_0_lvf = 0.1203;                   % LV ED pressure param (mmHg)
                P_0_rvf = 0.1203;                   % RV ED pressure param (mmHg)
                SVFact = 1.0000; % 1.3333           % Stress blood vol factor (uls)
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
            BestParams2Load = ['BCHF' num2str(RunPtNum) 'Red4_OptimBestParams.txt'];
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
        
    end   
            

%% **********************************************************************************
%  Opt/Sim Params of        R E D U C E D   S M I T H   M O D E L   S I M / O P T
% ***********************************************************************************

    % Simulation parameters
    NumBeats_SS = SimOptSpecs_Vect(1);                      % Num beats 2 stdy state
    NumBeats_ResPlot = SimOptSpecs_Vect(2);                 % Num beats to plot
    % Optimization parameters
    if (HandTune_Flag == 0)
        fmc_MaxFunEvals = SimOptSpecs_Vect(3);              % Max num fmc evals
        fmc_MaxIter = SimOptSpecs_Vect(4);                  % Max num fmc iters
        if (Optim_Flag == 2)
            ga_PopSize = SimOptSpecs_Vect(5);               % GA population size
            ga_MaxStallGens = SimOptSpecs_Vect(6);          % GA num stall gens to
        end                                                 %  terminate optim
    end

    % Pick the parameters we want to adjust with numbers shown below:
    
    % 1 --> E_es_lvf        LV free wall elast (mmHg/mL)        LEFT VENTRICLE
    % 2 --> P_0_lvf         LV ED pressure param (mmHg)
    % 3 --> lambda_lvf      LV ED pressure param (1/mL)
    % 4 --> E_es_rvf        RV free wall elast (mmHg/mL)        RIGHT VENTRICLE 
    % 5 --> P_0_rvf         RV ED pressure param (mmHg)
    % 6 --> lambda_rvf      RV ED pressure param (1/mL)
    % 7 --> E_es_pa         Pulm arterial elastance (mmHg/mL)   PULMONARY 
    % 8 --> E_es_pu         Pulm venous elastance (mmHg/mL)     VASCULATURE
    % 9 --> R_pul           Pulm vasc resist (mmHg*s/mL)
    % 10 -> E_es_sa         Systemic arter elstnce (mmHg/mL)    SYSTEMIC 
    % 11 -> E_es_sv         Systemic venous elstnce (mmHg/mL)   VASCULATURE
    % 12 -> R_sys           Syst art resistance (mmHg*s/mL)
    % 13 -> R_mt            Mitral valve resist (mmHg*s/mL)     HEART VALVES
    % 14 -> R_av            Aortic valve resist (mmHg*s/mL)
    % 15 -> R_tc            Tricspd vlv resist (mmHg*s/mL)
    % 16 -> R_pv            Pulmon vlv resist (mmHg*s/mL)
    % 17 -> SVFact          Stressed blood vol factor (uls)     CIRC BLOOD VOLUME
    
    % Bounds on all possible parameters
    LowBp_All(1) = log(0.100);      % LEFT          % LV free wall elast, E_es_lvf
    UpBp_All(1) = log(10.000);      %  VENTRICLE
    LowBp_All(2) = log(0.01);                       % LV ED pressure param, P_0_lvf
    UpBp_All(2) = log(5.00);
    LowBp_All(3) = log(0.005);                      % LV ED press param, lambda_lvf
    UpBp_All(3) = log(0.100);
    LowBp_All(4) = log(0.050);     % RIGHT          % RV free wall elast, E_es_rvf
    UpBp_All(4) = log(5.000);      %  VENTRICLE
    LowBp_All(5) = log(0.01);                       % RV ED pressure param, P_0_rvf
    UpBp_All(5) = log(5.00);
    LowBp_All(6) = log(0.005);                      % RV ED press param, lambda_rvf
    UpBp_All(6) = log(0.100);
    LowBp_All(7) = log(0.050);     % PULMONARY      % Pulm artery elastance, E_es_pa
    UpBp_All(7) = log(5.000);      %  VASCULATURE
    LowBp_All(8) = log(0.0005);                     % Pulm vein elastance, E_es_pu
    UpBp_All(8) = log(0.1000);
    LowBp_All(9) = log(0.005);                      % Pulm vasc resist, R_pul
    UpBp_All(9) = log(1.000);
    LowBp_All(10) = log(0.05);      % SYSTEMIC      % Aorta elastance, E_es_ao
    UpBp_All(10) = log(5.00);       %  VASCULATURE
    LowBp_All(11) = log(0.0001);                    % Vena cava elastance, E_es_vc
    UpBp_All(11) = log(0.1000);
    LowBp_All(12) = log(0.05);                      % Syst art resistance, R_sys
    UpBp_All(12) = log(15);
    LowBp_All(13) = log(0.005);     % HEART         % Mitral valve resist, R_mt
    UpBp_All(13) = log(0.500);      %  VALVES
    LowBp_All(14) = log(0.005);                     % Aortic valve resist, R_av
    UpBp_All(14) = log(0.500);
    LowBp_All(15) = log(0.005);                     % Tricspd vlv resist, R_tc
    UpBp_All(15) = log(0.500);
    LowBp_All(16) = log(0.0005);                    % Pulmon vlv resist, R_pv
    UpBp_All(16) = log(0.2500);
    LowBp_All(17) = log(0.1);       % CIRC          % Stress blood vol fact, SVFact
    UpBp_All(17) = log(3);          %  BLOOD VOL                          
    
    if (HandTune_Flag == 0 && NomParam_Flag == 0)
        
        % Create the adjustable parameter vector and parameter name strings 
        %  to pass to the objective function to so adjustable parameters can
        %  be written into CVParam_Struct on each optimization evaluation
        % Assign the adjustable parameter vector according to input file
        Num_AdjParams = size(AdjParam_Vect,1) - ...     % Num parameters
        sum(isnan(AdjParam_Vect));
        AdjParams = AdjParam_Vect(1:Num_AdjParams);     % Adj params up to first NaN
        
        p = zeros(1,Num_AdjParams);
        AdjParam_Strngs = cell(Num_AdjParams,1);
        LowBp = zeros(Num_AdjParams,1);
        UpBp = zeros(Num_AdjParams,1);
        for i = 1:Num_AdjParams
            ParamNum = AdjParams(i);
            p(i) = log(CVParam_Values{ParamNum+1});
            AdjParam_Strngs(i,1) = CVParam_Fields(ParamNum+1);
            LowBp(i) = LowBp_All(ParamNum);
            UpBp(i) = UpBp_All(ParamNum);
        end
        
    else
        
        AdjParam_Strngs = {};
        
    end

    % Building the simulation parameter structure
    SimOptParam_Values = {AdjParam_Strngs NumBeats_SS NumBeats_ResPlot};
    SimOptParam_Fields = {'AdjParam_Strngs' 'NumBeats_SS' 'NumBeats_ResPlot'};
    SimOptParam_Struct = cell2struct(SimOptParam_Values, ...
        SimOptParam_Fields,2);
    
    % Now put all structures into a single structure to pass
    if (NormParam_Flag == 1)
        AllStruct_Values = {FlagData_Struct PatData_Struct ...
            CVParam_Struct SimOptParam_Struct};
        AllStruct_Fields = {'FlagData_Struct' 'PatData_Struct' ...
            'CVParam_Struct' 'SimOptParam_Struct'};
        AllStruct_Struct = cell2struct(AllStruct_Values, ...
        AllStruct_Fields,2);
    else
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
    end
        
        
%% **********************************************************************************
%  Optimization of          R E D U C E D   S M I T H   M O D E L   S I M / O P T
% ***********************************************************************************
    
    if (HandTune_Flag == 0 && OptimBest_Flag == 0 && ...
            NormParam_Flag == 0 && NomParam_Flag == 0)
        
        % Setting fmincon options for fmincon optimization
        %  alone or the fmincon optimization after the 
        %  genetic algorithm optimization
            fmc_Opts = optimset('fmincon');             % Grab default settings
            fmc_Opts.Algorithm = 'active-set';          % Set algorithm type
            fmc_Opts.Display = 'iter';                  % Set display content
            fmc_Opts.TolFun = 1e-10;                    % Set tolerance on function
            fmc_Opts.MaxFunEvals = ... %5000;           % Set max num of fun evals
                fmc_MaxFunEvals; 
            fmc_Opts.MaxIter = fmc_MaxIter; %100;       % Set max num of iterations
            if (Parallel_Flag == 1)
                fmc_Opts.UseParallel = 1;
            end
            
        % PERFORMING OPTIMIZATION
        if (Optim_Flag == 1)
            
            % fmincon optimization call
            p_Optim = fmincon(@PatSpcHFRed4_ObjFun,p,[],[],[],[], ...
                LowBp,UpBp,[],fmc_Opts,AllStruct_Struct);
            
        elseif (Optim_Flag == 2 && Parallel_Flag == 0)
            
            % Genetic Algorithm Opimization - Serial computation
            ga_Opts = optimoptions('ga', ...
                'PopulationSize',ga_PopSize, ...            % 250
                'Display','iter', ...
                'MaxStallGenerations',ga_maxStallGens);     % 10
            Num_AdjParams = size(p,2);
            % Set objective function handle and 
            %  identify adjustable parameter vector, p
            PatSpecHFRed4_ObjFun_Hndl = ...
                @(p) PatSpecHFRed4_ObjFun(p,AllStruct_Struct);
            % Optimization call
            p_OptimGA = ga(PatSpecHFRed4_ObjFun_Hndl,Num_AdjParams,[],[],[],[], ...
                LowBp,UpBp,[],ga_Opts);
            p_Optim = fmincon(@PatSpecHFRed4_ObjFun,p_OptimGA,[],[],[],[], ...
                LowBp,UpBp,[],fmc_Opts,AllStruct_Struct);
            
        else
            
            % Genetic Algorithm Opimization - Parallel computation
            if (max(size(gcp)) == 0)                    % Parallel pool needed
                parpool                                 % Create the parallel pool
            end
            spmd
                warning('off','all');
            end  
            ga_Opts = optimoptions('ga', ...
                'PopulationSize',ga_PopSize, ...            % 250
                'Display','iter', ...
                'MaxStallGenerations',ga_MaxStallGens, ...  % 10
                'UseParallel',true);
            Num_AdjParams = size(p,2);
            % Set objective function handle and 
            %  identify adjustable parameter vector, p
            PatSpecHFRed4_ObjFun_Hndl = ...
                @(p) PatSpecHFRed4_ObjFun(p,AllStruct_Struct);
            % Optimization call
            [p_OptimGA,Res_OptimGA,~,~,FinPopGA,FinScoreGA] = ...
                ga(PatSpecHFRed4_ObjFun_Hndl,Num_AdjParams, ...
                [],[],[],[],LowBp,UpBp,[],ga_Opts);
            p_Optim = fmincon(@PatSpecHFRed4_ObjFun,p_OptimGA,[], ...
                [],[],[],LowBp,UpBp,[],fmc_Opts,AllStruct_Struct);
            
            % Overwrite initial adjustable parameter values with optimized values
            Num_AdjParams = size(p,2);
            CVParam_FieldNames = fieldnames(CVParam_Struct);
            Num_AllParams = size(CVParam_FieldNames,1);
            for i = 1:Num_AdjParams
                for j = 1:Num_AllParams
                    if (strcmp(CVParam_FieldNames{j},AdjParam_Strngs{i}))
                    CVParam_Struct.(AdjParam_Strngs{i}) = exp(p_Optim(i));
                    end
                end
            end
            
        end
            
    end
    
    
%% **********************************************************************************
%  Run Simulation of        R E D U C E D   S M I T H   M O D E L   S I M / O P T
% ***********************************************************************************
           
    if (NormParam_Flag == 1)                          % Normal parameter simulation
        
        % RUN SIMULATION OF REDUCED MODEL WITH NORMAL OR NOMINAL PARAMETERS 
        % Calculate total blood volume based on height, weight and sex. 
        %  This expression is from Nadler et al. Surgery 51:224,1962.
        if (Sex == 'M')
            TotBV = ((0.3669 * (Hgt/100)^3) + ...
                (0.03219 * BW) + 0.6041) * 1000;
        else
            TotBV = ((0.3561 * (Hgt/100)^3) + ...
                (0.03308 * BW) + 0.1833) * 1000;
        end
        % The original Smith model only circulated a portion of the blood so aortic
        %  pressure dynamics are not lumped into a general arterila systemic 
        %  compartment. Assuming they were simulating a typical 5000 mL total blood 
        %  volume they included only 1500 mL (or 30%) in the circulating volume
        %  therefore we will multiply our calculated TotBV value by 0.3 to yield 
        %  circulating blood volume. To account for extra recruited volume in heart
        %  disease the 30% circulating blood volume can be altered by changing SVFact
        SVFact = 1.00;                              % Smith Value of SVFact
        CircBV = SVFact * 0.30 * TotBV;
        % Setting state variable initial conditions 
        %  Note that initial volume division is scaled as a fraction
        %  of circulating blood volume calculated earlier and the 
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
        % Build driver function parameter structure
        HR = HR_Norm;                                 % Heart rate (beats/min)
        period = 60/HR;                               % Period of heart beat (s)
        B = HR;                                       % Elstnce fctn param (1/s^2)
        C = period/2;                                 % Elastance fctn param (s)
        DriverP_Values = {HR period B C};
        DriverP_Fields = {'HR' 'period' 'B' 'C'};
        DriverP_Struct = cell2struct(DriverP_Values, ...
            DriverP_Fields,2);
        % Calculating timespans to reach steady state and then for simulation
        TSpan_SS = [0 NumBeats_SS * period];
        % Solve over the time span with ode15s
        [T_Out_Norm,X_Out_Norm] = ode15s(@dXdT_SmithRed4, ...
            TSpan_SS,X0,[],DriverP_Struct,CVParam_Struct);
        % Calculation intermediate pressures to plot
        Num_TOut_Norm = size(T_Out_Norm,1); % Number of time points
        P_LV_Norm = zeros(Num_TOut_Norm,1); % Preallocating matrices
        P_RV_Norm = zeros(Num_TOut_Norm,1);
        P_SA_Norm = zeros(Num_TOut_Norm,1);
        P_SV_Norm = zeros(Num_TOut_Norm,1);
        P_PA_Norm = zeros(Num_TOut_Norm,1);
        P_PU_Norm = zeros(Num_TOut_Norm,1);
        % Running Smith model and loading intermediates for plotting
        for i = 1:Num_TOut_Norm
            VarOut = dXdT_SmithRed4(T_Out_Norm(i), ...
                X_Out_Norm(i,:),DriverP_Struct, ...
                CVParam_Struct,1);
            P_LV_Norm(i) = VarOut(1);
            P_RV_Norm(i) = VarOut(2);
            P_SA_Norm(i) = VarOut(3);
            P_SV_Norm(i) = VarOut(4);
            P_PA_Norm(i) = VarOut(5);
            P_PU_Norm(i) = VarOut(6);
        end
        
    else
        
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
        % The original Smith model only circulated a portion of the blood so aortic
        %  pressure dynamics are not lumped into a general arterila systemic 
        %  compartment. Assuming they were simulating a typical 5000 mL total blood 
        %  volume they included only 1500 mL (or 30%) in the circulating volume
        %  therefore we will multiply our calculated TotBV value by 0.3 to yield 
        %  circulating blood volume. To account for extra recruited volume in heart
        %  disease the 30% circulating blood volume can be altered by changing SVFact
        SVFact = CVParam_Struct.SVFact;
        CircBV = SVFact * 0.30 * TotBV;
        % Setting state variable initial conditions 
        %  Note that initial volume division is scaled as a fraction
        %  of circulating blood volume calculated earlier and the 
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
            
            plot(T_Out_RHC,X_Out_RHC(:,1))
            
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
            
            % GETTING THE RESIDUAL AT THE OPTIM PARAMETER VALUES
            if (HandTune_Flag == 0 && OptimBest_Flag == 0 && ...
                    NormParam_Flag == 0 && NomParam_Flag == 0)
                Res_Optim = PatSpecHFRed4_ObjFun(p_Optim,AllStruct_Struct);
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
        
            % GETTING THE RESIDUAL AT THE OPTIM PARAMETER VALUES
            if (HandTune_Flag == 0 && OptimBest_Flag == 0 && NormParam_Flag == 0)
                Res_Optim = PatSpecHFRed4_ObjFun(p_Optim,AllStruct_Struct);
            end
            
        end
        
    end
    

%% **********************************************************************************
%  Plot Sim or Opt of   S M I T H   C A R D I O V A S C   M O D E L   S I M / O P T
% ***********************************************************************************

    % PLOT NORMAL SIMULATION OR OPTIM/HAND TUNE WITH RHC/ECHO OR RHC ONLY DATA
    %  If we are just simulating the Smith et al. result then no data is 
    %  plotted and we run one simulation. If we have optimization results or 
    %  are performing a hand tune we plot out the simulation with a RHC simulation
    %  and an Echo simulation when both datasets are present and only a RHC
    %  simulation if only RHC data is present.
    
    if (NormParam_Flag == 1)
        
        % NORMAL CARDIOVASCULAR FUNCTION
        % Finding the index and the time associated with the last
        %  few beats to plot which were used for the residual calculation
        NumBeat_Start = NumBeats_SS - NumBeats_ResPlot;     % Starting beat
        tStart_Ind = find(T_Out_Norm >= ...                 % Grab the start index
            (NumBeat_Start * period),1,'first');
        tStart = T_Out_Norm(tStart_Ind);                    % Calc start time
        % Saving the the simulations of the last few beats to plot
        T_Out = T_Out_Norm(tStart_Ind:end) - tStart;        % Setting start to 0
        P_RVSim = P_RV_Norm(tStart_Ind:end);                % P_RV to plot
        P_SASim = P_SA_Norm(tStart_Ind:end);                % P_SA to plot
        P_PASim = P_PA_Norm(tStart_Ind:end);                % P_PA to plot
        P_PUSim = P_PU_Norm(tStart_Ind:end);                % P_PU to plot
        V_LVSim = X_Out_Norm(tStart_Ind:end,1);             % V_LV to plot
        V_RVSim = X_Out_Norm(tStart_Ind:end,2);             % V_RV to plot
        P_LVSim = P_LV_Norm(tStart_Ind:end);                % P_LV to plot
        P_SVSim = P_SV_Norm(tStart_Ind:end);                % P_SV to plot
        % Capturing min and max left ventricular volume of the 
        %  last few beats to calculate cardiac output
        V_LVES_Norm = 1000;
        V_LVED_Norm = 0;
        for i = tStart_Ind:Num_TOut_Norm
            V_LVES_Norm = ...
                min(V_LVES_Norm,X_Out_Norm(i,1));
            V_LVED_Norm = ...
                max(V_LVED_Norm,X_Out_Norm(i,1));
        end
        CO_Norm = ((V_LVED_Norm - ...
            V_LVES_Norm) * HR_Norm) / 1000;
        CO_Sim = CO_Norm;
        % Now storing everything into a structure to send to the
        %  plotting function SixPanel_Figure
        SimPlot_Values = {T_Out P_RVSim P_SASim P_PASim ...
            P_PUSim V_LVSim V_RVSim P_LVSim P_SVSim CO_Sim};
        SimPlot_Fields = {'T_Out' 'P_RVSim' 'P_SASim' 'P_PASim' ...
            'P_PUSim' 'V_LVSim' 'V_RVSim' 'P_LVSim' 'P_SVSim' 'CO_Sim'};
        SimPlot_Struct = cell2struct(SimPlot_Values, ...
            SimPlot_Fields,2);
        % Plotting function call
        SixPanel_Figure(AllStruct_Struct,SimPlot_Struct)
        
        % Saving normal cardiovascular variables so it can be used
        %  later to plot against the patient specific results
        if (SaveNorm_Flag == 1)
            % Saving the values under another name for later
            T_SimNorm = T_Out;                              % Normal CV time 
            P_RVSimNorm = P_RVSim;                          % Normal CV P_RV 
            P_SASimNorm = P_SASim;                          % Normal CV P_SA
            P_PASimNorm = P_PASim;                          % Normal CV P_PA
            P_PUSimNorm = P_PUSim;                          % Normal CV P_PU
            P_LVSimNorm = P_LVSim;                          % Normal CV P_LV
            P_SVSimNorm = P_SVSim;                          % Normal CV P_SV
            V_LVSimNorm = V_LVSim;                          % Normal CV V_LV
            V_RVSimNorm = V_RVSim;                          % Normal CV V_RV
            CO_SimNorm = CO_Sim;                            % Normal CV CO
            % Save in local folder as NormVars.mat
            SaveFile1 = 'NormVars.mat';
            save(SaveFile1, 'T_SimNorm', 'P_RVSimNorm', ...
                'P_SASimNorm', 'P_PASimNorm', 'P_PUSimNorm', ...
                'P_LVSimNorm', 'P_SVSimNorm','V_LVSimNorm', ...
                'V_RVSimNorm','CO_SimNorm')
        end
        
    else
        
        if (RHCEcho_Flag == 1)
            
            % PATIENT SPECIFIC MODEL RESULTS - ONLY RHC + TTE
            % Finding the index and the time associated with the last
            %  few beats to plot which were used for the residual calculation
            NumBeat_Start = NumBeats_SS - NumBeats_ResPlot; % Starting beat
%             T_Out_RHC
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
            
            % CHECK TO MAKE SURE WE HAVE REACHED STEADY STATE
            % GET VOLUMES TO PLOT
            V_RVED_Echo = max(V_RVEchoSim);
            V_RVES_Echo = min(V_RVEchoSim); 
            CO_RVEchoSim = ((V_RVED_Echo - V_RVES_Echo) * HR_Echo) / 1000;
            VTotEcho = X_Out_Echo(tStartEcho_Ind:end,1) + ...
                X_Out_Echo(tStartEcho_Ind:end,2) + ...
                X_Out_Echo(tStartEcho_Ind:end,3) + ...
                X_Out_Echo(tStartEcho_Ind:end,4) + ...
                X_Out_Echo(tStartEcho_Ind:end,5) + ...
                X_Out_Echo(tStartEcho_Ind:end,6); 
            VTotRHC = X_Out_RHC(tStartRHC_Ind:end,1) + ...
                X_Out_RHC(tStartRHC_Ind:end,2) + ...
                X_Out_RHC(tStartRHC_Ind:end,3) + ...
                X_Out_RHC(tStartRHC_Ind:end,4) + ...
                X_Out_RHC(tStartRHC_Ind:end,5) + ...
                X_Out_RHC(tStartRHC_Ind:end,6); 
            
            % PLOT ALL IN ONE FIGURE
            ScrSize = get(0,'ScreenSize');              % Getting screen size
            VolCheck_Fig = figure('Position', ...       % Positioning the figure
                [ScrSize(3)/50 ScrSize(4)/50 ...        %  on the screen
                ScrSize(3)/2 ScrSize(4)/2]); 
            plot(T_OutRHC,VTotRHC,'-b','LineWidth',3)
            hold on
            plot(T_OutEcho,VTotEcho,':k','LineWidth',3)
            plot(T_OutEcho,X_Out_Echo(tStartEcho_Ind:end,1),'--r','LineWidth',3)
            plot(T_OutEcho,X_Out_Echo(tStartEcho_Ind:end,2),'--g','LineWidth',3)
            plot(T_OutEcho,X_Out_Echo(tStartEcho_Ind:end,3),'--b','LineWidth',3)
            plot(T_OutEcho,X_Out_Echo(tStartEcho_Ind:end,4),':r','LineWidth',3)
            plot(T_OutEcho,X_Out_Echo(tStartEcho_Ind:end,5),':g','LineWidth',3)
            plot(T_OutEcho,X_Out_Echo(tStartEcho_Ind:end,6),':b','LineWidth',3)
            hold off
            legend('V_{Tot,RHC}','V_{Tot,Echo}','V_{LV,Echo}', ...
                'V_{RV,Echo}', 'V_{PA,Echo}', 'V_{PU,Echo}', ...
                'V_{AO,Echo}', 'V_{VC,Echo}')
            xlabel('Time, s','FontSize',14,'FontWeight','bold')
            ylabel('Volume, mL','FontSize',14,'FontWeight','bold')
            
            
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
        
    end
    
    % Saving optimized parameter values, the adjusted parameters
    %  and the residual with the optimized parameters
    if (HandTune_Flag == 0 && NormParam_Flag == 0 && ...
            OptimBest_Flag == 0 && NomParam_Flag == 0)
        Expp_Optim = exp(p_Optim);
        SaveFile2 = 'OptimParams.mat';
        save(SaveFile2,'Expp_Optim','AdjParam_Strngs','Res_Optim')
    end

    toc