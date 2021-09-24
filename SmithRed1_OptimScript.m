% ***********************************************************************************
%          S M I T H   C A R D I O V A S C U L A R   S Y S T E M S   M O D E L
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
%   ventricular-ventricular interaction has been removed.
%
%   Model originally created on     17 January 2016
%   Model last modfied on           10 January 2019

%   Developed by        Brian Carlson
%                       Physiological Systems Dynamics Laboratory
%                       Department of Molecular and Integrative Physiology
%                       University of Michigan
%
% ***********************************************************************************
%  Start of             S M I T H   C A R D I O V A S C   M O D E L   S I M / O P T
% ***********************************************************************************

%% **********************************************************************************
%  Flag Data for        S M I T H   C A R D I O V A S C   M O D E L   S I M / O P T
% ***********************************************************************************

    tic
    warning('off','all');
    
    DataStruct = importdata('BCHF18_InputData.txt');
    
    FlagData_Vect = DataStruct.data(1:8,1);
    PatData_Vect = DataStruct.data(1:4,2);
    RHCData_Vect = DataStruct.data(1:10,3);
    if (size(DataStruct.data,2) == 4)
        EchoData_Vect = DataStruct.data(1:4,4);
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
    SmithParam_Flag = FlagData_Vect(6);         % Run Smith params = 1, not = 0
    SaveSmith_Flag = FlagData_Vect(7);          % Save Smith run = 1, not = 0
    Comp2Smith_Flag = FlagData_Vect(8);         % Compare to Smith = 1, not = 0   
    
    % Save all flag data so it can be passed to functions
    FlagData_Values = {HandTune_Flag OptimBest_Flag Optim_Flag ...
        Parallel_Flag RHCEcho_Flag SmithParam_Flag ...
        Comp2Smith_Flag SaveSmith_Flag};
    FlagData_Fields = {'HandTune_Flag' 'OptimBest_Flag' ...
        'Optim_Flag' 'Parallel_Flag' 'RHCEcho_Flag' ...
        'SmithParam_Flag' 'Comp2Smith_Flag' 'SaveSmith_Flag'};
    FlagData_Struct = cell2struct(FlagData_Values, ...
                FlagData_Fields,2);
            
            
%% **********************************************************************************
%  Patient Data for     S M I T H   C A R D I O V A S C   M O D E L   S I M / O P T
% ***********************************************************************************
            
    if (SmithParam_Flag == 0)
            
        % We are running patient data simulation or optimization so get the 
        % general data of study patient number, weight, height and gender
        DID_Num = ...                               % Deidentified subject number
            ['BCHF' num2str(PatData_Vect(1))];                     
        BW = PatData_Vect(2);                       % Body weight (kg)
        Hgt = PatData_Vect(3);                      % Height (cm)
        if (PatData_Vect(4) == 1)                   % Gender (M or F)
            Gender = 'M';  
        else
            Gender = 'F';
        end
        % Save all patient data into a structure to be passed to functions
        PatData_Values = {DID_Num BW Hgt Gender};
        PatData_Fields = {'DID_Num' 'BW' 'Hgt' 'Gender'};
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
       	CO_Thermo = RHCData_Vect(10);               % Cardiac outpt thermodil (L/min)
        % Save all RHC data into a structure to be passed to functions
        RHCData_Values = {P_RVsyst P_RVdiast P_PAsyst P_PAdiast ....
            P_PCWave P_AOsyst P_AOdiast HR_RHC CO_Fick CO_Thermo};
        RHCData_Fields = {'P_RVsyst' 'P_RVdiast' 'P_PAsyst' 'P_PAdiast' ....
            'P_PCWave' 'P_AOsyst' 'P_AOdiast' 'HR_RHC' 'CO_Fick' 'CO_Thermo'};
        RHCData_Struct = cell2struct(RHCData_Values, ...
            RHCData_Fields,2);
        
        % Get the echo data if present in the input file
        if (size(DataStruct.data,2) == 4)
            V_LVsyst = EchoData_Vect(1);            % Systolic LV volume (mL)
            V_LVdiast = EchoData_Vect(2);           % Diastolic LV volume (mL)
            HR_Echo	= EchoData_Vect(3);             % Average heart rate (beats/min)
            CO_EchoD = EchoData_Vect(4);            % Cardiac output Echo-Dop (L/min)
            % Save all Echo data into a structure to be passed to functions
            EchoData_Values = {V_LVsyst V_LVdiast HR_Echo CO_EchoD};
            EchoData_Fields = {'V_LVsyst' 'V_LVdiast' 'HR_Echo' 'CO_EchoD'};
            EchoData_Struct = cell2struct(EchoData_Values, ...
                EchoData_Fields,2);
        end
        	
    else
        	
        DID_Num = 'Normal';
        HR_Smith = 80;                              % Average heart rate (beats/min)
        BW = 72.3;                                  % Body weight (kg, ~160 lbs)
        Hgt = 178;                                  % Height (cm, ~5'10")
        Gender = 'M';                               % Assuming male since TotBV ~ 5L
        
        PatData_Values = {DID_Num HR_Smith BW Hgt Gender};
        PatData_Fields = {'DID_Num' 'HR_Smith' 'BW' 'Hgt' 'Gender'};
        PatData_Struct = cell2struct(PatData_Values, ...
            PatData_Fields,2);
        
    end
    

%% **********************************************************************************
%  Fixed Params of      S M I T H   C A R D I O V A S C   M O D E L   S I M / O P T
% ***********************************************************************************
    
    if (SmithParam_Flag == 1)
         
        % SMITH ET AL. MODEL PARAMETERS
        % Elastance function driver parameters
        A = 1;                                      % Elastance function param (uls)
        % Pericardium calculation parameters
        P_0_pcd = 0.5003;                           % Unstress pericard press (mmHg)
        lambda_pcd = 0.03;                          % Pericard press exp term (1/mL)
        V_0_pcd = 200;                              % Pericard zero P volume (mL)
        P_th = -4;                                  % Thoracic pressure (mmHg)
        % Left ventricle free wall parameters
        E_es_lvf = 2.8798;                          % LV free wall elast (mmHg/mL) 
        V_d_lvf = 0;                                % LV ES zero P volume (mL)        
        P_0_lvf = 0.1203;                           % LV ED pressure param (mmHg)
        lambda_lvf = 0.033;                         % LV ED pressure param (1/mL)
        V_0_lvf = 0;                                % LV ED pressure param (mL)    
        % Right ventricle free wall parameters
        E_es_rvf = 0.585;                           % RV free wall elast (mmHg/mL) 
        V_d_rvf = 0;                                % RV ES zero P volume (mL)
        P_0_rvf = 0.2157;                           % RV ED pressure param (mmHg)
        lambda_rvf = 0.023;                         % RV ED pressure param (1/mL)
        V_0_rvf = 0;                                % RV ED pressure param (mL)
        % Pulmonary artery and vein parameters
        E_es_pa = 0.369;                            % Pulm artery elastance (mmHg/mL)
        V_d_pa = 0;                                 % Pulm artery zero P volume (mL)
        E_es_pu = 0.0073;                           % Pulm vein elastance (mmHg/mL)
        V_d_pu = 0;                                 % Pulm vein zero P volume (mL)
        R_pul = 0.1552;                             % Pulm vasc resist (mmHg*s/mL)
        % Aortic and vena cava parameters
        E_es_ao = 0.6913;                           % Aorta elastance (mmHg/mL)
        V_d_ao = 0;                                 % Aorta zero P volume (mL)
        E_es_vc = 0.0059;                           % Vena cava elastance (mmHg/mL)
        V_d_vc = 0;                                 % Vena cava zero P volume (mL)
        R_sys = 1.0889;                             % Syst art resistance (mmHg*s/mL)
        % Heart valve paramenters
        R_mt = 0.0158;                              % Mitral valve resist (mmHg*s/mL)
        L_mt = 7.6968e-5;                           % Mitrl vlve inert (mmHg*s^2/mL)
        R_av = 0.018;                               % Aortic valve resist (mmHg*s/mL)
        L_av = 1.2189e-4;                           % Aortic vlve inert (mmHg*s^2/mL)
        R_tc = 0.0237;                              % Tricspd vlv resist (mmHg*s/mL)
        L_tc = 8.0093e-5;                           % Tricspd vlv inert (mmHg*s^2/mL)
        R_pv = 0.0055;                              % Pulmon vlv resist (mmHg*s/mL)
        L_pv = 1.4868e-4;                           % Pulmon vlv inert (mmHg*s^2/mL)
        % Septum free wall parameters
        E_es_spt = 48.754;                          % Septum FW elstnce (mmHg/mL)
        V_d_spt = 2;                                % Septum zero P volume (mL)
        P_0_spt = 1.1101;                           % Septum ED pressure param (mmHg)
        lambda_spt = 0.435;                         % Septum ED pressure param (1/mL)
        V_0_spt = 2;                                % Septum ED pressure param (mL)
        % Heart failure param
        SVFact = 1.00;                              % Stress blood vol factor (uls)
        
        % Save all parameters in a sturcture to pass
        CVParam_Values = {A P_0_pcd lambda_pcd V_0_pcd P_th E_es_lvf ...
            V_d_lvf P_0_lvf lambda_lvf V_0_lvf E_es_rvf V_d_rvf P_0_rvf ...
            lambda_rvf V_0_rvf E_es_pa V_d_pa E_es_pu V_d_pu R_pul E_es_ao ...
            V_d_ao E_es_vc V_d_vc R_sys R_mt L_mt R_av L_av R_tc L_tc R_pv L_pv ...
            E_es_spt V_d_spt P_0_spt lambda_spt V_0_spt SVFact};
        CVParam_Fields = {'A' 'P_0_pcd' 'lambda_pcd' 'V_0_pcd' 'P_th' ...
            'E_es_lvf' 'V_d_lvf' 'P_0_lvf' 'lambda_lvf' 'V_0_lvf' 'E_es_rvf' ...
            'V_d_rvf' 'P_0_rvf' 'lambda_rvf' 'V_0_rvf' 'E_es_pa' 'V_d_pa' ...
            'E_es_pu' 'V_d_pu' 'R_pul' 'E_es_ao' 'V_d_ao' 'E_es_vc' 'V_d_vc' ...
            'R_sys' 'R_mt' 'L_mt' 'R_av' 'L_av' 'R_tc' 'L_tc' 'R_pv' 'L_pv' ...
            'E_es_spt' 'V_d_spt' 'P_0_spt' 'lambda_spt' 'V_0_spt' 'SVFact'};
        CVParam_Struct = cell2struct(CVParam_Values,CVParam_Fields,2);
        
    else
        
        if (OptimBest_Flag == 0)
            
            % MODEL PARAMETERS TO ADJUST
            % Elastance function driver parameters
            A = 1;                           %FIX   % Elastance function param (uls)
            % Pericardium calculation parameters
            P_0_pcd = 0.5003;                %FIX   % Unstress pericard press (mmHg)
            lambda_pcd = 0.03;               %FIX   % Pericard press exp term (1/mL)
            V_0_pcd = 200*2.00;              %FIX   % Pericard zero P volume (mL)
            P_th = 0; %-4;              %ADJ        % Thoracic pressure (mmHg)
            % Left ventricle free wall parameters
            E_es_lvf = 2.8798*0.575;    %ADJ        % LV free wall elast (mmHg/mL) 
            V_d_lvf = 0;                %ADJ        % LV ES zero P volume (mL)        
            P_0_lvf = 0.1203;           %ADJ        % LV ED pressure param (mmHg)
            lambda_lvf = 0.033*0.5575;  %ADJ        % LV ED pressure param (1/mL)
            V_0_lvf = 0;                %ADJ        % LV ED pressure param (mL)    
            % Right ventricle free wall parameters
            E_es_rvf = 0.585*0.75;      %ADJ        % RV free wall elast (mmHg/mL) 
            V_d_rvf = 0;                %ADJ        % RV ES zero P volume (mL)
            P_0_rvf = 0.2157;           %ADJ        % RV ED pressure param (mmHg)
            lambda_rvf = 0.023*0.55;    %ADJ        % RV ED pressure param (1/mL)
            V_0_rvf = 0;                %ADJ        % RV ED pressure param (mL)
            % Pulmonary artery and vein parameters
            E_es_pa = 0.369*0.70;       %ADJ        % Pulm artery elastance (mmHg/mL)
            V_d_pa = 0;                      %FIX   % Pulm artery zero P volume (mL)
            E_es_pu = 0.0073*75.00;     %ADJ        % Pulm vein elastance (mmHg/mL)
            V_d_pu = 0;                      %FIX   % Pulm vein zero P volume (mL)
            R_pul = 0.1552*1.190;       %ADJ        % Pulm vasc resist (mmHg*s/mL)
            % Aortic and vena cava parameters
            E_es_ao = 0.6913*1.65;      %ADJ        % Aorta elastance (mmHg/mL)
            V_d_ao = 0;                      %FIX   % Aorta zero P volume (mL)
            E_es_vc = 0.0059;           %ADJ        % Vena cava elastance (mmHg/mL)
            V_d_vc = 0;                      %FIX   % Vena cava zero P volume (mL)
            R_sys = 1.0889*0.9375;      %ADJ        % Syst art resistance (mmHg*s/mL)
            % Heart valve paramenters
            R_mt = 0.0158;                   %FIX   % Mitral valve resist (mmHg*s/mL)
            L_mt = 7.6968e-5;                %FIX   % Mitrl vlve inert (mmHg*s^2/mL)
            R_av = 0.018;                    %FIX   % Aortic valve resist (mmHg*s/mL)
            L_av = 1.2189e-4;                %FIX   % Aortic vlve inert (mmHg*s^2/mL)
            R_tc = 0.0237;                   %FIX   % Tricspd vlv resist (mmHg*s/mL)
            L_tc = 8.0093e-5;                %FIX   % Tricspd vlv inert (mmHg*s^2/mL)
            R_pv = 0.0055*2.00;         %ADJ        % Pulmon vlv resist (mmHg*s/mL)
            L_pv = 1.4868e-4;                %FIX   % Pulmon vlv inert (mmHg*s^2/mL)
            % Heart failure param
            SVFact = 1.00;              %ADJ        % Stress blood vol factor (uls)
            
            % Save all parameters in a sturcture to pass
            CVParam_Values = {A P_0_pcd lambda_pcd V_0_pcd P_th ...
                E_es_lvf V_d_lvf P_0_lvf lambda_lvf V_0_lvf E_es_rvf ...
                V_d_rvf P_0_rvf lambda_rvf V_0_rvf E_es_pa V_d_pa ...
                E_es_pu V_d_pu R_pul E_es_ao V_d_ao E_es_vc V_d_vc ...
                R_sys R_mt L_mt R_av L_av R_tc L_tc R_pv L_pv SVFact};
            CVParam_Fields = {'A' 'P_0_pcd' 'lambda_pcd' 'V_0_pcd' ...
                'P_th' 'E_es_lvf' 'V_d_lvf' 'P_0_lvf' 'lambda_lvf' ...
                'V_0_lvf' 'E_es_rvf' 'V_d_rvf' 'P_0_rvf' 'lambda_rvf' ...
                'V_0_rvf' 'E_es_pa' 'V_d_pa' 'E_es_pu' 'V_d_pu' ...
                'R_pul' 'E_es_ao' 'V_d_ao' 'E_es_vc' 'V_d_vc' 'R_sys' ...
                'R_mt' 'L_mt' 'R_av' 'L_av' 'R_tc' 'L_tc' 'R_pv' ...
                'L_pv' 'SVFact'};
            CVParam_Struct = cell2struct(CVParam_Values,CVParam_Fields,2);
        
        else
        
            % LOAD THE OPTIM PARAMS FROM AN EARLIER RUN
             % Elastance function driver parameters
            A = 1;                                  % Elastance function param (uls)
            % Pericardium calculation parameters
            P_0_pcd = 0.5003;                       % Unstress pericard press (mmHg)
            lambda_pcd = 0.03;                      % Pericard press exp term (1/mL)
            V_0_pcd = 400;                          % Pericard zero P volume (mL)
            P_th = -0.553369764822130;              % Thoracic pressure (mmHg)
            % Left ventricle free wall parameters
            E_es_lvf = exp(-0.448877316482953);     % LV free wall elast (mmHg/mL) 
            V_d_lvf = exp(-12.937640357765114);     % LV ES zero P volume (mL)        
            P_0_lvf = exp(-2.955665215375613);      % LV ED pressure param (mmHg)
            lambda_lvf = exp(-3.133231196988598);   % LV ED pressure param (1/mL)
            V_0_lvf = exp(-9.074023887724229);      % LV ED pressure param (mL)    
            % Right ventricle free wall parameters
            E_es_rvf = exp(-2.055413438087283);     % RV free wall elast (mmHg/mL) 
            V_d_rvf = exp(-1.904585759433639);      % RV ES zero P volume (mL)
            P_0_rvf = exp(-4.317117340318973);      % RV ED pressure param (mmHg)
            lambda_rvf = exp(-3.930459013835557);   % RV ED pressure param (1/mL)
            V_0_rvf = exp(4.506423412520583);       % RV ED pressure param (mL)
            % Pulmonary artery and vein parameters
            E_es_pa = exp(0.578271494753880);       % Pulm artery elastance (mmHg/mL)
            V_d_pa = 0;                             % Pulm artery zero P volume (mL)
            E_es_pu = exp(-2.368138759177784);      % Pulm vein elastance (mmHg/mL)
            V_d_pu = 0;                             % Pulm vein zero P volume (mL)
            R_pul = exp(-2.795688940630368);        % Pulm vasc resist (mmHg*s/mL)
            % Aortic and vena cava parameters
            E_es_ao = exp(0.577007037053961);       % Aorta elastance (mmHg/mL)
            V_d_ao = 0;                             % Aorta zero P volume (mL)
            E_es_vc = exp(-3.910877999228612);      % Vena cava elastance (mmHg/mL)
            V_d_vc = 0;                             % Vena cava zero P volume (mL)
            R_sys = exp(0.534077537829140);         % Syst art resistance (mmHg*s/mL)
            % Heart valve paramenters
            R_mt = 0.0158;                          % Mitral valve resist (mmHg*s/mL)
            L_mt = 7.6968e-5;                       % Mitrl vlve inert (mmHg*s^2/mL)
            R_av = 0.018;                           % Aortic valve resist (mmHg*s/mL)
            L_av = 1.2189e-4;                       % Aortic vlve inert (mmHg*s^2/mL)
            R_tc = 0.0237;                          % Tricspd vlv resist (mmHg*s/mL)
            L_tc = 8.0093e-5;                       % Tricspd vlv inert (mmHg*s^2/mL)
            R_pv = exp(-3.357757662731381);         % Pulmon vlv resist (mmHg*s/mL)
            L_pv = 1.4868e-4;                       % Pulmon vlv inert (mmHg*s^2/mL)
            % Heart failure param
            SVFact = exp(0.535240661598561);        % Stress blood vol factor (uls)
            
            % Save all parameters in a sturcture to pass
            CVParam_Values = {A P_0_pcd lambda_pcd V_0_pcd P_th ...
                E_es_lvf V_d_lvf P_0_lvf lambda_lvf V_0_lvf E_es_rvf ...
                V_d_rvf P_0_rvf lambda_rvf V_0_rvf E_es_pa V_d_pa ...
                E_es_pu V_d_pu R_pul E_es_ao V_d_ao E_es_vc V_d_vc ...
                R_sys R_mt L_mt R_av L_av R_tc L_tc R_pv L_pv SVFact};
            CVParam_Fields = {'A' 'P_0_pcd' 'lambda_pcd' 'V_0_pcd' ...
                'P_th' 'E_es_lvf' 'V_d_lvf' 'P_0_lvf' 'lambda_lvf' ...
                'V_0_lvf' 'E_es_rvf' 'V_d_rvf' 'P_0_rvf' 'lambda_rvf' ...
                'V_0_rvf' 'E_es_pa' 'V_d_pa' 'E_es_pu' 'V_d_pu' ...
                'R_pul' 'E_es_ao' 'V_d_ao' 'E_es_vc' 'V_d_vc' 'R_sys' ...
                'R_mt' 'L_mt' 'R_av' 'L_av' 'R_tc' 'L_tc' 'R_pv' ...
                'L_pv' 'SVFact'};
            CVParam_Struct = cell2struct(CVParam_Values,CVParam_Fields,2);
            
        end
        
    end

        
%% **********************************************************************************
%  Opt/Sim Params of    S M I T H   C A R D I O V A S C   M O D E L   S I M / O P T
% ***********************************************************************************
        
    % Pick the parameters we want to adjust with numbers shown below:
    
    % 1 --> P_0_pcd         Unstress pericard press (mmHg)      PERICARDIUM and
    % 2 --> lambda_pcd      Pericard press exp term (1/mL)      THORACIC CHAMBER
    % 3 --> V_0_pcd         Pericard zero P volume (mL)
    % 4 --> P_th            Thoracic pressure (mmHg)
    % 5 --> E_es_lvf        LV free wall elast (mmHg/mL)        LEFT VENTRICLE
    % 6 --> V_d_lvf         LV ES zero P volume (mL)        
    % 7 --> P_0_lvf         LV ED pressure param (mmHg)
    % 8 --> lambda_lvf      LV ED pressure param (1/mL)
    % 9 --> V_0_lvf         LV ED pressure param (mL)    
    % 10 -> E_es_rvf        RV free wall elast (mmHg/mL)        RIGHT VENTRICLE 
    % 11 -> V_d_rvf         RV ES zero P volume (mL)
    % 12 -> P_0_rvf         RV ED pressure param (mmHg)
    % 13 -> lambda_rvf      RV ED pressure param (1/mL)
    % 14 -> V_0_rvf         RV ED pressure param (mL)
    % 15 -> E_es_pa         Pulm artery elastance (mmHg/mL)     PULMONARY 
    % 16 -> V_d_pa          Pulm artery zero P volume (mL)      VASCULATURE
    % 17 -> E_es_pu         Pulm vein elastance (mmHg/mL)
    % 18 -> V_d_pu          Pulm vein zero P volume (mL)
    % 19 -> R_pul           Pulm vasc resist (mmHg*s/mL)
    % 20 -> E_es_ao         Aorta elastance (mmHg/mL)           SYSTEMIC 
    % 21 -> V_d_ao          Aorta zero P volume (mL)            VASCULATURE
    % 22 -> E_es_vc         Vena cava elastance (mmHg/mL)
    % 23 -> V_d_vc          Vena cava zero P volume (mL)
    % 24 -> R_sys           Syst art resistance (mmHg*s/mL)
    % 25 -> R_mt            Mitral valve resist (mmHg*s/mL)     HEART VALVES
    % 26 -> L_mt            Mitrl vlve inert (mmHg*s^2/mL)
    % 27 -> R_av            Aortic valve resist (mmHg*s/mL)
    % 28 -> L_av            Aortic vlve inert (mmHg*s^2/mL)
    % 29 -> R_tc            Tricspd vlv resist (mmHg*s/mL)
    % 30 -> L_tc            Tricspd vlv inert (mmHg*s^2/mL)
    % 31 -> R_pv            Pulmon vlv resist (mmHg*s/mL)
    % 32 -> L_pv            Pulmon vlv inert (mmHg*s^2/mL)
    % 33 -> SVFact          Stressed blood vol factor (uls)     CIRC BLOOD VOLUME
    
    % SELECT ADJUSTABLE PARAMETERS HERE
    AdjParams = [4;5;6;7;8;9;10;11;12; ...
        13;14;15;17;19;20;22;24;31;33];
    Num_AdjParams = size(AdjParams,1);
    
    % Bounds on all possible parameters
    LowBp_All(1) = log(0.01);       % PERICARD &    % Unstr pericard P, P_0_pcd
    UpBp_All(1) = log(5.00);        %  THORACIC
    LowBp_All(2) = log(0.005);                      % Pericard press exp, lambda_pcd
    UpBp_All(2) = log(0.100);
    LowBp_All(3) = log(100);                        % Pericard zero P vol, V_0_pcd
    UpBp_All(3) = log(700);
    LowBp_All(4) = -5;                              % Thoracic pressure, P_th
    UpBp_All(4) = 5;
    LowBp_All(5) = log(0.100);      % LEFT          % LV free wall elast, E_es_lvf
    UpBp_All(5) = log(10.000);      %  VENTRICLE
    LowBp_All(6) = log(0.000001);                   % LV ES zero P volume, V_d_lvf
    UpBp_All(6) = log(150);
    LowBp_All(7) = log(0.01);                       % LV ED pressure param, P_0_lvf
    UpBp_All(7) = log(5.00);
    LowBp_All(8) = log(0.005);                      % LV ED press param, lambda_lvf
    UpBp_All(8) = log(0.100);
    LowBp_All(9) = log(0.000001);                   % LV ED pressure param, V_0_lvf
    UpBp_All(9) = log(150);
    LowBp_All(10) = log(0.050);     % RIGHT         % RV free wall elast, E_es_rvf
    UpBp_All(10) = log(5.000);      %  VENTRICLE
    LowBp_All(11) = log(0.0001);                    % RV ES zero P volume, V_d_rvf
    UpBp_All(11) = log(150);
    LowBp_All(12) = log(0.01);                      % RV ED pressure param, P_0_rvf
    UpBp_All(12) = log(5.00);
    LowBp_All(13) = log(0.005);                     % RV ED press param, lambda_rvf
    UpBp_All(13) = log(0.100);
    LowBp_All(14) = log(0.000001);                   % RV ED pressure param, V_0_rvf
    UpBp_All(14) = log(150);
    LowBp_All(15) = log(0.050);     % PULMONARY     % Pulm artery elastance, E_es_pa
    UpBp_All(15) = log(5.000);      %  VASCULATURE
    LowBp_All(16) = log(0.000001);                  % Pulm artery zero P vol, V_d_pa
    UpBp_All(16) = log(10);
    LowBp_All(17) = log(0.0005);                    % Pulm vein elastance, E_es_pu
    UpBp_All(17) = log(0.1000);
    LowBp_All(18) = log(0.000001);                  % Pulm vein zero P vol, V_d_pu
    UpBp_All(18) = log(15);
    LowBp_All(19) = log(0.005);                     % Pulm vasc resist, R_pul
    UpBp_All(19) = log(1.000);
    LowBp_All(20) = log(0.05);      % SYSTEMIC      % Aorta elastance, E_es_ao
    UpBp_All(20) = log(5.00);       %  VASCULATURE
    LowBp_All(21) = log(0.000001);                  % Aorta zero P volume, V_d_ao
    UpBp_All(21) = log(20);
    LowBp_All(22) = log(0.0001);                    % Vena cava elastance, E_es_vc
    UpBp_All(22) = log(0.1000);
    LowBp_All(23) = log(0.000001);                  % Vena cava zero P vol, V_d_vc
    UpBp_All(23) = log(30);
    LowBp_All(24) = log(0.05);                      % Syst art resistance, R_sys
    UpBp_All(24) = log(15);
    LowBp_All(25) = log(0.005);     % HEART         % Mitral valve resist, R_mt
    UpBp_All(25) = log(0.500);      %  VALVES
    LowBp_All(26) = log(0.5e-5);                    % Mitrl vlve inert, L_mt
    UpBp_All(26) = log(150.0e-5);
    LowBp_All(27) = log(0.005);                     % Aortic valve resist, R_av
    UpBp_All(27) = log(0.500);
    LowBp_All(28) = log(0.05e-4);                   % Aortic valve inert, L_av
    UpBp_All(28) = log(100.00e-4);
    LowBp_All(29) = log(0.005);                     % Tricspd vlv resist, R_tc
    UpBp_All(29) = log(0.500);
    LowBp_All(30) = log(0.05e-5);                   % Tricspd vlv inert, L_tc
    UpBp_All(30) = log(150.00e-5);
    LowBp_All(31) = log(0.0005);                    % Pulmon vlv resist, R_pv
    UpBp_All(31) = log(0.2500);
    LowBp_All(32) = log(0.05e-4);                   % Pulmon vlv inert, L_pv
    UpBp_All(32) = log(100.00e-4);
    LowBp_All(33) = log(0.1);       % CIRC          % Stress blood vol fact, SVFact
    UpBp_All(33) = log(3);          %  BLOOD VOL                          
    
    if (HandTune_Flag == 0)
        
        % Create the adjustable parameter vector and parameter name strings 
        %  to pass to the objective function to so adjustable parameters can
        %  be written into CVParam_Struct on each optimization evaluation
        p = zeros(1,Num_AdjParams);
        AdjParam_Strngs = cell(Num_AdjParams,1);
        LowBp = zeros(Num_AdjParams,1);
        UpBp = zeros(Num_AdjParams,1);
        for i = 1:Num_AdjParams
            ParamNum = AdjParams(i);
            if (strcmp(CVParam_Fields(ParamNum+1),'P_th'))
                p(i) = CVParam_Values{ParamNum+1};
            else
                p(i) = log(CVParam_Values{ParamNum+1});
            end
            AdjParam_Strngs(i,1) = CVParam_Fields(ParamNum+1);
            LowBp(i) = LowBp_All(ParamNum);
            UpBp(i) = UpBp_All(ParamNum);
        end
        
    else
        
        AdjParam_Strngs = {};
        
    end

    % Number of beats to steady state and then for sim/optim
    NumBeats_SS = 20;
    NumBeats_Sim = 10;
    NumBeats_PlotRes = 5;
    
    % Building the simulation parameter structure
    SimOptParam_Values = {AdjParam_Strngs NumBeats_SS ...
        NumBeats_Sim NumBeats_PlotRes};
    SimOptParam_Fields = {'AdjParam_Strngs' 'NumBeats_SS' ...
        'NumBeats_Sim' 'NumBeats_PlotRes'};
    SimOptParam_Struct = cell2struct(SimOptParam_Values, ...
        SimOptParam_Fields,2);
    
    % Now put all structures into a single structure to pass
    if (SmithParam_Flag == 1)
        AllStruct_Values = {FlagData_Struct PatData_Struct ...
            CVParam_Struct SimOptParam_Struct};
        AllStruct_Fields = {'FlagData_Struct' 'PatData_Struct' ...
            'CVParam_Struct' 'SimOptParam_Struct'};
        AllStruct_Struct = cell2struct(AllStruct_Values, ...
        AllStruct_Fields,2);
    else
        if (size(DataStruct.data,2) == 4)
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
%  Optimization of      S M I T H   C A R D I O V A S C   M O D E L   S I M / O P T
% ***********************************************************************************
    
    if (HandTune_Flag == 0 && OptimBest_Flag == 0 && SmithParam_Flag == 0)
        
        % Setting fmincon options for fmincon optimization
        %  alone or the fmincon optimization after the 
        %  genetic algorithm optimization
            fmc_Opts = optimset('fmincon');             % Grab default settings
            fmc_Opts.Algorithm = 'active-set';          % Set algorithm type
            fmc_Opts.Display = 'iter';                  % Set display content
            fmc_Opts.TolFun = 1e-10;                    % Set tolerance on function
            fmc_Opts.MaxFunEvals = 150;   % 2000        % Set max num of fun evals
            fmc_Opts.MaxIter = 10;        % 100         % Set max num of iterations
            if (Parallel_Flag == 1)
                fmc_Opts.UseParallel = 1;
            end
            
        % PERFORMING OPTIMIZATION
        if (Optim_Flag == 1)
            
            % fmincon optimization call
            p_Optim = fmincon(@PatSpecHF_SmithRed1ObjFun,p,[],[], ...
                [],[],LowBp,UpBp,[],fmc_Opts,AllStruct_Struct);
            
        elseif (Optim_Flag == 2 && Parallel_Flag == 0)
            
            % Genetic Algorithm Opimization - Serial computation
            ga_Opts = optimoptions('ga', ...
                'PopulationSize',15, ...            % 250
                'Display','iter', ...
                'MaxStallGenerations',3, ...        % 10
                'OutputFcn', @GAGenStore);
            Num_AdjParams = size(p,2);
            % Set objective function handle and 
            %  identify adjustable parameter vector, p
            PatSpecHF_SmithRed1ObjFun_Hndl = ...
                @(p) PatSpecHF_SmithRed1ObjFun(p,AllStruct_Struct);
            % Optimization call
            p_OptimGA = ga(PatSpecHF_SmithRed1ObjFun_Hndl, ...
                Num_AdjParams,[],[],[],[],LowBp,UpBp,[],ga_Opts);
            p_Optim = fmincon(@PatSpecHF_SmithRed1ObjFun,p_OptimGA, ...
                [],[],[],[],LowBp,UpBp,[],fmc_Opts,AllStruct_Struct);
            
        else
            
            % Genetic Algorithm Opimization - Parallel computation
            if (max(size(gcp)) == 0)                    % Parallel pool needed
                parpool                                 % Create the parallel pool
            end
            spmd
                warning('off','all');
            end  
            ga_Opts = optimoptions('ga', ...
                'PopulationSize',15, ...            % 250
                'Display','iter', ...
                'MaxStallGenerations',3, ...        % 10
                'OutputFcn', @GAGenStore, ...
                'UseParallel',true);
            Num_AdjParams = size(p,2);
            % Set objective function handle and 
            %  identify adjustable parameter vector, p
            PatSpecHF_SmithRed1ObjFun_Hndl = ...
                @(p) PatSpecHF_SmithRed1ObjFun(p,AllStruct_Struct);
            % Optimization call
            p_OptimGA = ga(PatSpecHF_SmithRed1ObjFun_Hndl, ...
                Num_AdjParams,[],[],[],[],LowBp,UpBp,[],ga_Opts);
            p_Optim = fmincon(@PatSpecHF_SmithRed1ObjFun,p_OptimGA, ...
                [],[],[],[],LowBp,UpBp,[],fmc_Opts,AllStruct_Struct);
            
        end
        
        % Overwrite initial adjustable parameter values with optimized values
        Num_AdjParams = size(p,2);
        CVParam_FieldNames = fieldnames(CVParam_Struct);
        Num_AllParams = size(CVParam_FieldNames,1);
        for i = 1:Num_AdjParams
            for j = 1:Num_AllParams
                if (strcmp(CVParam_FieldNames{j},AdjParam_Strngs{i}))
                    if (strcmp(AdjParam_Strngs{i},'P_th'))
                        CVParam_Struct.(AdjParam_Strngs{i}) = p_Optim(i);
                    else
                        CVParam_Struct.(AdjParam_Strngs{i}) = exp(p_Optim(i));
                    end
                end
            end
        end
            
    end
    
    
    
    
%% **********************************************************************************
%  Run Simulation of    S M I T H   C A R D I O V A S C   M O D E L   S I M / O P T
% ***********************************************************************************
    
    if (SmithParam_Flag == 1)                       % Smith et al. simulation
        
        % RUN SMITH ET AL. SIMULATION
        % Calculate total blood volume based on height, weight and gender. 
        %  This expression is from Nadler et al. Surgery 51:224,1962.
        if (Gender == 'M')
            TotBV = ((0.3669 * (Hgt/100)^3) + ...
                (0.03219 * BW) + 0.6041) * 1000;
        else
            TotBV = ((0.3561 * (Hgt/100)^3) + ...
                (0.03308 * BW) + 0.1833) * 1000;
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
        % Put into vector to pass to ode15s
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
        % Build driver function parameter structure
        HR = HR_Smith;                              % RHC heart rate (beats/min)
        period = 60/HR_Smith;                       % Period of heart beat (s)
        B = HR_Smith;                               % Elastance fctn param (1/s^2)
        C = period/2;                               % Elastance fctn param (s)
        DriverP_Values = {HR period B C};
        DriverP_Fields = {'HR' 'period' 'B' 'C'};
        DriverP_Struct = cell2struct(DriverP_Values, ...
            DriverP_Fields,2);
        % Calculating the initial condition on the volume due to the septum 
        %  deflection which is an implicit function requiring the use of fsolve
        time = 0;
        FSolve_Opts = optimset('Diagnostics','off','Display','off');
        SeptZF_Hndl = @(V_spt) SeptZF(V_spt,V_lv0, ...
            V_rv0,time,DriverP_Struct,CVParam_Struct);
        [V_spt0,SeptZF_Val,ExitFlag,OutData] = ...
            fsolve(SeptZF_Hndl,-15, FSolve_Opts);
        % Put last initial condition into vector
        X0(11) = V_spt0;
        % Calculating timespans to reach steady state and then for simulation
        TSpan_SS = [0 NumBeats_SS * period];
        TSpan_Sim = [0 NumBeats_Sim * period];
        % Build mass matrix for DAE solver
        M = eye(11);                                % Put identity on diagonal
        M(11,11) = 0;                               % Set last express as a ZFun
        % Set ODE/DAE options and time span
        ODE_Opts = odeset('Mass',M); 
        % Solve over the time span with ode15s
        [T_Out_SmithSS,X_Out_SmithSS] = ode15s(@dXdT_Smith, ...
            TSpan_SS,X0,ODE_Opts,DriverP_Struct,CVParam_Struct);
        % Now solve over the simulation time span
        [T_Out_SmithSim,X_Out_SmithSim] = ode15s(@dXdT_Smith,TSpan_Sim, ...
            X_Out_SmithSS(end,:),ODE_Opts,DriverP_Struct,CVParam_Struct);
        % Calculation intermediate pressures to plot
        Num_TOut_SmithSim = size(T_Out_SmithSim,1); % Number of time points
        P_LV_SmithSim = zeros(Num_TOut_SmithSim,1); % Preallocating matrices
        P_RV_SmithSim = zeros(Num_TOut_SmithSim,1);
        P_AO_SmithSim = zeros(Num_TOut_SmithSim,1);
        P_VC_SmithSim = zeros(Num_TOut_SmithSim,1);
        P_PA_SmithSim = zeros(Num_TOut_SmithSim,1);
        P_PU_SmithSim = zeros(Num_TOut_SmithSim,1);
        V_LVES_SmithSim = 1000;
        V_LVED_SmithSim = 0;
        % Running Smith model and loading intermediates for plotting
        for i = 1:Num_TOut_SmithSim
            VarOut = dXdT_Smith(T_Out_SmithSim(i), ...
                X_Out_SmithSim(i,:),DriverP_Struct, ...
                CVParam_Struct,1);
            P_LV_SmithSim(i) = VarOut(1);
            P_RV_SmithSim(i) = VarOut(2);
            P_AO_SmithSim(i) = VarOut(3);
            P_VC_SmithSim(i) = VarOut(4);
            P_PA_SmithSim(i) = VarOut(5);
            P_PU_SmithSim(i) = VarOut(6);
            V_LVES_SmithSim = ...
                min(V_LVES_SmithSim,X_Out_SmithSim(i,1));
            V_LVED_SmithSim = ...
                max(V_LVED_SmithSim,X_Out_SmithSim(i,1));
        end
        CO_SmithSim = ((V_LVED_SmithSim - ...
            V_LVES_SmithSim) * HR_Smith) / 1000;
        
    else                                            % Optim or handtune
        
        % Calculate total blood volume based on height, weight and gender
        %  which is the same for both RHC and Echo simulations
        %  This expression is from Nadler et al. Surgery 51:224,1962.
        if (Gender == 'M')
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
        V_ao0 = (133.338/1500) * CircBV;
        V_vc0 = (329.780/1500) * CircBV;
        Q_mt0 = TotBV/60;
        Q_av0 = TotBV/60;
        Q_tc0 = TotBV/60;
        Q_pv0 = TotBV/60;
        % Put into vector to pass to ode15s
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
        
        if (RHCEcho_Flag == 1)                      % With RHC and Echo data
            
            % RUN RHC SIMULATION FIRST
            % Build driver function parameter structure
            HR = HR_RHC;                            % RHC heart rate (beats/min)
            period = 60/HR_RHC;                     % Period of heart beat (s)
            B = HR_RHC;                             % Elastance fctn param (1/s^2)
            C = period/2;                           % Elastance fctn param (s)
            DriverP_Values = {HR period B C};
            DriverP_Fields = {'HR' 'period' 'B' 'C'};
            DriverP_Struct = cell2struct(DriverP_Values, ...
                DriverP_Fields,2);
            %Calculating timespans to reach steady state and then for simulation
            TSpan_SS = [0 NumBeats_SS * period];
            TSpan_Sim = [0 NumBeats_Sim * period];
            % Solve over the steady state time span with ode15s
            [T_Out_RHCSS,X_Out_RHCSS] = ode15s(@dXdT_SmithRed1, ...
                TSpan_SS,X0,[],DriverP_Struct,CVParam_Struct);
            % Now solve over the simulation time span with ode15s
            [T_Out_RHCSim,X_Out_RHCSim] = ode15s(@dXdT_SmithRed1, ...
                TSpan_Sim,X_Out_RHCSS(end,:),[],DriverP_Struct,CVParam_Struct);
            
            % CALCULATING INTERMEDIATE PRESSURES TO PLOT
            Num_TOut_RHCSim = size(T_Out_RHCSim,1); % Number of time points
            P_LV_RHCSim = zeros(Num_TOut_RHCSim,1); % Preallocating matrices
            P_RV_RHCSim = zeros(Num_TOut_RHCSim,1);
            P_AO_RHCSim = zeros(Num_TOut_RHCSim,1);
            P_VC_RHCSim = zeros(Num_TOut_RHCSim,1);
            P_PA_RHCSim = zeros(Num_TOut_RHCSim,1);
            P_PU_RHCSim = zeros(Num_TOut_RHCSim,1);
            V_LVES_RHCSim = 1000;
            V_LVED_RHCSim = 0;
            % RUNNING MODEL TO GET INTERMEDIATES
            for i = 1:Num_TOut_RHCSim
                VarOut = dXdT_SmithRed1(T_Out_RHCSim(i), ...
                    X_Out_RHCSim(i,:),DriverP_Struct, ...
                    CVParam_Struct,1);
                P_LV_RHCSim(i) = VarOut(1);
                P_RV_RHCSim(i) = VarOut(2);
                P_AO_RHCSim(i) = VarOut(3);
                P_VC_RHCSim(i) = VarOut(4);
                P_PA_RHCSim(i) = VarOut(5);
                P_PU_RHCSim(i) = VarOut(6);
                V_LVES_RHCSim = ...
                    min(V_LVES_RHCSim,X_Out_RHCSim(i,1));
                V_LVED_RHCSim = ...
                    max(V_LVED_RHCSim,X_Out_RHCSim(i,1));
            end
            CO_RHCSim = ((V_LVED_RHCSim - ...
                V_LVES_RHCSim) * HR_RHC) / 1000;
            
            % RUN ECHO SIMULATION NEXT
            % Build driver function parameter structure
            HR = HR_Echo;                           % RHC heart rate (beats/min)
            period = 60/HR_Echo;                    % Period of heart beat (s)
            B = HR_Echo;                            % Elastance fctn param (1/s^2)
            C = period/2;                           % Elastance fctn param (s)
            DriverP_Values = {HR period B C};
            DriverP_Fields = {'HR' 'period' 'B' 'C'};
            DriverP_Struct = cell2struct(DriverP_Values, ...
                DriverP_Fields,2);
            % Solve over the steady state time span with ode15s
            [T_Out_EchoSS,X_Out_EchoSS] = ode15s(@dXdT_SmithRed1, ...
                TSpan_SS,X0,[],DriverP_Struct,CVParam_Struct);
            % Now solve over the simulation time span with ode15s
            [T_Out_EchoSim,X_Out_EchoSim] = ode15s(@dXdT_SmithRed1, ...
                TSpan_Sim,X_Out_EchoSS(end,:),[], ...
                DriverP_Struct,CVParam_Struct);
            
            % CALCULATING INTERMEDIATE PRESSURES TO PLOT
            Num_TOut_EchoSim = ...                  % Number of time points
                size(T_Out_EchoSim,1); 
            P_LV_EchoSim = ...                      % Preallocating matrices
                zeros(Num_TOut_EchoSim,1);   
            P_RV_EchoSim = ...
                zeros(Num_TOut_EchoSim,1);
            % RUNNING MODEL TO GET INTERMEDIATES
            for i = 1:Num_TOut_EchoSim
                VarOut = dXdT_SmithRed1(T_Out_EchoSim(i), ...
                    X_Out_EchoSim(i,:),DriverP_Struct, ...
                    CVParam_Struct,1);
                P_LV_EchoSim(i) = VarOut(1);
                P_RV_EchoSim(i) = VarOut(2);
            end
            
            % CALCULATING CO OF ECHO SIMULATION
            V_LVES_EchoSim = 1000;
            V_LVED_EchoSim = 0;
            % FINDING LV END DIASTOLIC AND END SYSTOLIC VOLUME
            for i = 1:Num_TOut_EchoSim
                V_LVES_EchoSim = ...
                    min(V_LVES_EchoSim,X_Out_EchoSim(i,1));
                V_LVED_EchoSim = ...
                    max(V_LVED_EchoSim,X_Out_EchoSim(i,1));
            end
            CO_EchoSim = ((V_LVED_EchoSim - ...
                V_LVES_EchoSim) * HR_Echo) / 1000;
            
            % GETTING THE RESIDUAL AT THE OPTIM PARAMETER VALUES
            if (HandTune_Flag == 0 && OptimBest_Flag == 0 && SmithParam_Flag == 0)
                Res_Optim = PatSpecHF_SmithRed1ObjFun(p_Optim,AllStruct_Struct);
            end
            
        else                                        % RHC data only
            
            % RUN RHC SIMULATION ONLY
            % Build driver function parameter structure
            HR = HR_RHC;                            % RHC heart rate (beats/min)
            period = 60/HR_RHC;                     % Period of heart beat (s)
            B = HR_RHC;                             % Elastance fctn param (1/s^2)
            C = period/2;                           % Elastance fctn param (s)
            DriverP_Values = {HR period B C};
            DriverP_Fields = {'HR' 'period' 'B' 'C'};
            DriverP_Struct = cell2struct(DriverP_Values, ...
                DriverP_Fields,2);
            %Calculating timespans to reach steady state and then for simulation
            TSpan_SS = [0 NumBeats_SS * period];
            TSpan_Sim = [0 NumBeats_Sim * period];
            % Solve over the steady state time span with ode15s
            [T_Out_RHCSS,X_Out_RHCSS] = ode15s(@dXdT_SmithRed1, ...
                TSpan_SS,X0,ODE_Opts,DriverP_Struct,CVParam_Struct);
            % Now solve over the simulation time span with ode15s
            [T_Out_RHCSim,X_Out_RHCSim] = ode15s(@dXdT_SmithRed1,TSpan_Sim, ...
                X_Out_RHCSS(end,:),ODE_Opts,DriverP_Struct,CVParam_Struct);
            
            % CALCULATING INTERMEDIATE PRESSURES TO PLOT
            Num_TOut_RHCSim = size(T_Out_RHCSim,1); % Number of time points
            P_LV_RHCSim = zeros(Num_TOut_RHCSim,1); % Preallocating matrices
            P_RV_RHCSim = zeros(Num_TOut_RHCSim,1);
            P_AO_RHCSim = zeros(Num_TOut_RHCSim,1);
            P_VC_RHCSim = zeros(Num_TOut_RHCSim,1);
            P_PA_RHCSim = zeros(Num_TOut_RHCSim,1);
            P_PU_RHCSim = zeros(Num_TOut_RHCSim,1);
            V_LVES_RHCSim = 1000;
            V_LVED_RHCSim = 0;
            % RUNNING MODEL TO GET INTERMEDIATES
            for i = 1:Num_TOut_RHCSim
                VarOut = dXdT_SmithRed1(T_Out_RHCSim(i), ...
                    X_Out_RHCSim(i,:),DriverP_Struct, ...
                    CVParam_Struct,1);
                P_LV_RHCSim(i) = VarOut(1);
                P_RV_RHCSim(i) = VarOut(2);
                P_AO_RHCSim(i) = VarOut(3);
                P_VC_RHCSim(i) = VarOut(4);
                P_PA_RHCSim(i) = VarOut(5);
                P_PU_RHCSim(i) = VarOut(6);
                V_LVES_RHCSim = ...
                    min(V_LVES_RHCSim,X_Out_RHCSim(i,1));
                V_LVED_RHCSim = ...
                    max(V_LVED_RHCSim,X_Out_RHCSim(i,1));
            end
            CO_RHCSim = ((V_LVED_RHCSim - ...
                V_LVES_RHCSim) * HR_RHC) / 1000;
            
            % GETTING THE RESIDUAL AT THE OPTIM PARAMETER VALUES
            if (HandTune_Flag == 0 && OptimBest_Flag == 0 && SmithParam_Flag == 0)
                Res_Optim = PatSpecHF_SmithRed1ObjFun(p_Optim,AllStruct_Struct);
            end
            
        end
        
    end
    
    
%% **********************************************************************************
%  Plot Sim or Opt of   S M I T H   C A R D I O V A S C   M O D E L   S I M / O P T
% ***********************************************************************************

    % PLOT SMITH SIMULATION OR OPTIM/HAND TUNE WITH RHC/ECHO OR RHC ONLY DATA
    %  If we are just simulating the Smith et al. result then no data is 
    %  plotted and we run one simulation. If we have optimization results or 
    %  are performing a hand tune we plot out the simulation with a RHC simulation
    %  and an Echo simulation when both datasets are present and only a RHC
    %  simulation if only RHC data is present.
    
    if (SmithParam_Flag == 1)
        
        T_Out = T_Out_SmithSim;
        P_RVSim = P_RV_SmithSim;
        P_AOSim = P_AO_SmithSim;
        P_PASim = P_PA_SmithSim;
        P_PUSim = P_PU_SmithSim;
        CO_Sim = CO_SmithSim;
        V_LVSim = X_Out_SmithSim(:,1);
        V_RVSim = X_Out_SmithSim(:,2);
        P_LVSim = P_LV_SmithSim;
        P_VCSim = P_VC_SmithSim;
        
        SimPlot_Values = {T_Out P_RVSim P_AOSim P_PASim ...
            P_PUSim CO_Sim V_LVSim V_RVSim P_LVSim P_VCSim};
        SimPlot_Fields = {'T_Out' 'P_RVSim' 'P_AOSim' 'P_PASim' ...
            'P_PUSim' 'CO_Sim' 'V_LVSim' 'V_RVSim' 'P_LVSim' 'P_VCSim'};
        SimPlot_Struct = cell2struct(SimPlot_Values, ...
            SimPlot_Fields,2);
     
        SixPanel_Figure(AllStruct_Struct,SimPlot_Struct)

        if (SaveSmith_Flag == 1)
            
            % Specifying the number of beats at the end of the RHC simulation
            %  that we want to calculate the residual from and then finding 
            %  the portion of the simulation to extract values from
            TimeStart_SmithPlot = ...
                (NumBeats_Sim - NumBeats_PlotRes) * period;
            tIndStart_SmithPlot = ...
                find((T_Out >= TimeStart_SmithPlot),1,'first');
            % Now store the last set of beats for plotting/saving
            T_SimSmith = ...
                T_Out(tIndStart_SmithPlot:end) - ...
                T_Out(tIndStart_SmithPlot);
            P_RVSimSmith = P_RVSim(tIndStart_SmithPlot:end);
            P_AOSimSmith = P_AOSim(tIndStart_SmithPlot:end);
            P_PASimSmith = P_PASim(tIndStart_SmithPlot:end);
            P_PUSimSmith = P_PUSim(tIndStart_SmithPlot:end);
            CO_SimSmith = CO_Sim;
            P_LVSimSmith = P_LVSim(tIndStart_SmithPlot:end);
            P_VCSimSmith = P_VCSim(tIndStart_SmithPlot:end);
            V_LVSimSmith = V_LVSim(tIndStart_SmithPlot:end);
            V_RVSimSmith = V_RVSim(tIndStart_SmithPlot:end);

            SaveFile1 = 'SmithVars.mat';
            save(SaveFile1, 'T_SimSmith', 'P_RVSimSmith', ...
                'P_AOSimSmith', 'P_PASimSmith', 'P_PUSimSmith', ...
                'CO_SimSmith', 'P_LVSimSmith', 'P_VCSimSmith',...
                'V_LVSimSmith', 'V_RVSimSmith')

        end
        
    else
        
        if (RHCEcho_Flag == 1)
            
            % Specifying the number of beats at the end of the RHC simulation
            %  that we want to calculate the residual from and then finding 
            %  the portion of the simulation to extract values from
            TimeStart_RHCPlot = ...
                (NumBeats_Sim - NumBeats_PlotRes) * period;
            tIndStart_RHCPlot = ...
                find((T_Out_RHCSim >= TimeStart_RHCPlot),1,'first');
            TimeStart_EchoPlot = ...
                (NumBeats_Sim - NumBeats_PlotRes) * period;
            tIndStart_EchoPlot = ...
                find((T_Out_EchoSim >= TimeStart_EchoPlot),1,'first');
            % Now store the last set of beats for plotting/saving
            T_OutRHC = ...
                T_Out_RHCSim(tIndStart_RHCPlot:end) - ...
                T_Out_RHCSim(tIndStart_RHCPlot);
            T_OutEcho = ...
                T_Out_EchoSim(tIndStart_EchoPlot:end) - ...
                T_Out_EchoSim(tIndStart_EchoPlot);
            P_RVRHCSim = P_RV_RHCSim(tIndStart_RHCPlot:end);
            P_RVEchoSim = P_RV_EchoSim(tIndStart_EchoPlot:end);
            P_AOSim = P_AO_RHCSim(tIndStart_RHCPlot:end);
            P_PASim = P_PA_RHCSim(tIndStart_RHCPlot:end);
            P_PUSim = P_PU_RHCSim(tIndStart_RHCPlot:end);
            V_LVSim = X_Out_EchoSim(tIndStart_EchoPlot:end,1);
            V_RVSim = X_Out_EchoSim(tIndStart_EchoPlot:end,2);
            P_LVRHCSim = P_LV_RHCSim(tIndStart_RHCPlot:end);
            P_LVEchoSim = P_LV_EchoSim(tIndStart_EchoPlot:end);
            P_VCSim = P_VC_RHCSim(tIndStart_RHCPlot:end);
        
            SimPlot_Values = {T_OutRHC T_OutEcho P_RVRHCSim ...
                P_RVEchoSim P_AOSim P_PASim P_PUSim ...
                CO_RHCSim CO_EchoSim V_LVSim V_RVSim ...
                P_LVRHCSim P_LVEchoSim P_VCSim};
            SimPlot_Fields = {'T_OutRHC' 'T_OutEcho' 'P_RVRHCSim' ...
                'P_RVEchoSim' 'P_AOSim' 'P_PASim' 'P_PUSim' ...
                'CO_RHCSim' 'CO_EchoSim' 'V_LVSim' 'V_RVSim' ...
                'P_LVRHCSim' 'P_LVEchoSim' 'P_VCSim'};
            SimPlot_Struct = cell2struct(SimPlot_Values, ...
                SimPlot_Fields,2);
     
            SixPanel_Figure(AllStruct_Struct,SimPlot_Struct)
            
        else
            
            % Specifying the number of beats at the end of the RHC simulation
            %  that we want to calculate the residual from and then finding 
            %  the portion of the simulation to extract values from
            TimeStart_RHCPlot = ...
                (NumBeats_Sim - NumBeats_PlotRes) * period;
            tIndStart_RHCPlot = ...
                find((T_Out_RHCSim >= TimeStart_RHCPlot),1,'first');
            % Now store the last set of beats for plotting/saving
            T_Out = ...
                T_Out_RHCSim(tIndStart_RHCPlot:end) - ...
                T_Out_RHCSim(tIndStart_RHCPlot);
            P_RVSim = P_RV_RHCSim(tIndStart_RHCPlot:end);
            P_AOSim = P_AO_RHCSim(tIndStart_RHCPlot:end);
            P_PASim = P_PA_RHCSim(tIndStart_RHCPlot:end);
            P_PUSim = P_PU_RHCSim(tIndStart_RHCPlot:end);
            V_LVSim = X_Out_RHCSim(tIndStart_RHCPlot:end,1);
            V_RVSim = X_Out_RHCSim(tIndStart_RHCPlot:end,2);
            P_LVSim = P_LV_RHCSim(tIndStart_RHCPlot:end);
            P_VCSim = P_VC_RHCSim(tIndStart_RHCPlot:end);
        
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
        
    if (HandTune_Flag == 0 && SmithParam_Flag == 0 && OptimBest_Flag == 0)
        Expp_Optim = exp(p_Optim);
        SaveFile2 = 'OptimParams.mat';
        save(SaveFile2,'Expp_Optim','AdjParam_Strngs','Res_Optim')
    end
    
    toc