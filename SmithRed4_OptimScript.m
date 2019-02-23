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
%   ventricular-ventricular interaction, valve inertances, all dead space/zero 
%   pressure volumes, pericardium and thoracic chambers have been removed.
%
%   Model originally created on     17   January 2016
%   Model last modfied on           20 September 2018

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
    
    % First we determine what we want to do with this code: optimize to data,
    %  compare to Smith model simulation, hand fit, run Smith model params and
    %  save output, or plot out a previous optimized fit. Specific detail on
    %  flag values for each option are given above in the script header
    RHCEcho_Flag = 1;                           % RHC + Echo = 1, RHC only = 2
    HandTune_Flag = 1;                          % 1 = Hnd tne or sim, 0 = Opt & sim
    Optim_Flag = 0;                             % 1 = fmincon, 2 = Genetic alg
    Parallel_Flag = 0;                          % 1 = Parallel, 0 = Serial optim
    OptimBest_Flag = 0;                         % 1 = Optim p Run 0 = Opt or Other
    SmithParam_Flag = 0;                        % Run Smith params = 1, not = 0
    SaveSmith_Flag = 0;                         % Save Smith run = 1, not = 0
    Comp2Smith_Flag = 1;                        % Compare to Smith = 1, not = 0   
    
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
            
    if (RHCEcho_Flag == 1 && SmithParam_Flag == 0)
            
        DID_Num = 'BCHF1';                      % Deidentified subject number
        P_RVsyst = 46;                          % Systolic RV pressure (mmHg)     *
        P_RVdiast = 10;                         % Diastolic RV pressure (mmHg)    *
        P_PAsyst = 44;                          % Systolic pulm art press (mmHg)  *
       	P_PAdiast = 20;                         % Diast pulm art press (mmHg)     *
       	P_PCWave = 19;                          % Average pulm wedge press (mmHg) * 
       	P_AOsyst = 166;                         % Systolic aortic press (mmHg)    *
       	P_AOdiast = 65;                         % Diastolic aortic press (mmHg)   *
       	Ave_HR = 56;                            % Average heart rate (beats/min)
       	CO_Fick = 6.13;                         % Cardiac output Fick (L/min)     *
       	CO_Thermo = 6.10;                       % Cardiac outpt thermodil (L/min) *
       	V_LVsyst = 99.0;                        % Systolic LV volume (mL)         *
        V_LVdiast = 210.0;                      % Diastolic LV volume (mL)        *
        BW = 107.4;                             % Body weight (kg)
        Hgt = 160;                              % Height (cm)
        Gender = 'F';                           % Gender (M or F)
            
        PatData_Values = {RHCEcho_Flag DID_Num P_RVsyst ...
            P_RVdiast P_PAsyst P_PAdiast P_PCWave P_AOsyst ...
            P_AOdiast Ave_HR CO_Fick CO_Thermo V_LVsyst ...
            V_LVdiast BW Hgt Gender};
        PatData_Fields = {'RHCEcho_Flag' 'DID_Num' 'P_RVsyst' ...
            'P_RVdiast' 'P_PAsyst' 'P_PAdiast' 'P_PCWave' ...
            'P_AOsyst' 'P_AOdiast' 'Ave_HR' 'CO_Fick' 'CO_Thermo' ...
            'V_LVsyst' 'V_LVdiast' 'BW' 'Hgt' 'Gender'};
        PatData_Struct = cell2struct(PatData_Values, ...
            PatData_Fields,2);
        	
    elseif (RHCEcho_Flag == 0 && SmithParam_Flag == 0)
        	
        DID_Num = 'BCHF1';                      % Deidentified subject number
		P_RVsyst = 46;                          % Systolic RV pressure (mmHg)     *
        P_RVdiast = 10;                         % Diastolic RV pressure (mmHg)    *
        P_PAsyst = 44;                          % Systolic pulm art press (mmHg)  *
        P_PAdiast = 20;                         % Diast pulm art press (mmHg)     *
        P_PCWave = 19;                          % Average pulm wedge press (mmHg) * 
        P_AOsyst = 166;                         % Systolic aortic press (mmHg)    *
        P_AOdiast = 65;                         % Diastolic aortic press (mmHg)   *
        Ave_HR = 56;                            % Average heart rate (beats/min)
        CO_Fick = 6.13;                         % Cardiac output Fick (L/min)     *
        CO_Thermo = 6.10;                       % Cardiac outpt thermodil (L/min) *
        BW = 107.4;                             % Body weight (kg)
        Hgt = 160;                              % Height (cm)
        Gender = 'F';                           % Gender (M or F)
        
        PatData_Values = {RHCEcho_Flag DID_Num P_RVsyst ...
            P_RVdiast P_PAsyst P_PAdiast P_PCWave P_AOsyst ...
            P_AOdiast Ave_HR CO_Fick CO_Thermo BW Hgt Gender};
        PatData_Fields = {'RHCEcho_Flag' 'DID_Num' 'P_RVsyst' ...
            'P_RVdiast' 'P_PAsyst' 'P_PAdiast' 'P_PCWave' ...
            'P_AOsyst' 'P_AOdiast' 'Ave_HR' 'CO_Fick' ...
            'CO_Thermo' 'BW' 'Hgt' 'Gender'};
        PatData_Struct = cell2struct(PatData_Values, ...
            PatData_Fields,2);
        
    else
        
        DID_Num = 'Normal';
        Ave_HR = 80;                            % Average heart rate (beats/min)
        BW = 72.3;                              % Body weight (kg, ~160 lbs)
        Hgt = 178;                              % Height (cm, ~5'10")
        Gender = 'M';                           % Assuming male since TotBV ~ 5L
        
        PatData_Values = {DID_Num Ave_HR BW Hgt Gender};
        PatData_Fields = {'DID_Num' 'Ave_HR' 'BW' 'Hgt' 'Gender'};
        PatData_Struct = cell2struct(PatData_Values, ...
            PatData_Fields,2);
        
    end
    

%% **********************************************************************************
%  Fixed Params of      S M I T H   C A R D I O V A S C   M O D E L   S I M / O P T
% ***********************************************************************************
    
    % Elastance function driver parameters
    period = 60/Ave_HR;                             % Period of heart beat (s)
    A = 1;                                          % Elastance function param (uls)
    B = Ave_HR;                                     % Elastance fctn param (1/s^2)
    C = period/2;                                   % Elastance fctn param (s)
    
    if (SmithParam_Flag == 1)
         
        % SMITH ET AL. MODEL PARAMETERS
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
        CVParam_Values = {period A B C P_0_pcd lambda_pcd V_0_pcd P_th E_es_lvf ...
            V_d_lvf P_0_lvf lambda_lvf V_0_lvf E_es_rvf V_d_rvf P_0_rvf ...
            lambda_rvf V_0_rvf E_es_pa V_d_pa E_es_pu V_d_pu R_pul E_es_ao ...
            V_d_ao E_es_vc V_d_vc R_sys R_mt L_mt R_av L_av R_tc L_tc R_pv L_pv ...
            E_es_spt V_d_spt P_0_spt lambda_spt V_0_spt SVFact};
        CVParam_Fields = {'period' 'A' 'B' 'C' 'P_0_pcd' 'lambda_pcd' 'V_0_pcd' ...
            'P_th' 'E_es_lvf' 'V_d_lvf' 'P_0_lvf' 'lambda_lvf' 'V_0_lvf' ...
            'E_es_rvf' 'V_d_rvf' 'P_0_rvf' 'lambda_rvf' 'V_0_rvf' 'E_es_pa' ...
            'V_d_pa' 'E_es_pu' 'V_d_pu' 'R_pul' 'E_es_ao' 'V_d_ao' 'E_es_vc' ...
            'V_d_vc' 'R_sys' 'R_mt' 'L_mt' 'R_av' 'L_av' 'R_tc' 'L_tc' 'R_pv' ...
            'L_pv' 'E_es_spt' 'V_d_spt' 'P_0_spt' 'lambda_spt' 'V_0_spt' 'SVFact'};
        CVParam_Struct = cell2struct(CVParam_Values,CVParam_Fields,2);
        
    else
        
        if (OptimBest_Flag == 0)
            
            % MODEL PARAMETERS TO ADJUST
            % Left ventricle free wall parameters
            E_es_lvf = 2.8798 * 0.580;              % LV free wall elast (mmHg/mL) 
            P_0_lvf = 0.1203;                       % LV ED pressure param (mmHg)
            lambda_lvf = 0.033 * 0.465;              % LV ED pressure param (1/mL)
            % Right ventricle free wall parameters
            E_es_rvf = 0.585 * 0.75;                % RV free wall elast (mmHg/mL) 
            P_0_rvf = 0.2157;                       % RV ED pressure param (mmHg)
            lambda_rvf = 0.023 * 0.525;              % RV ED pressure param (1/mL)
            % Pulmonary artery and vein parameters
            E_es_pa = 0.369 * 0.70;                 % Pulm artery elastance (mmHg/mL)
            E_es_pu = 0.0073 * 75.00;               % Pulm vein elastance (mmHg/mL)
            R_pul = 0.1552 * 1.35;                  % Pulm vasc resist (mmHg*s/mL)
            % Aortic and vena cava parameters
            E_es_ao = 0.6913 * 1.675;               % Aorta elastance (mmHg/mL)
            E_es_vc = 0.0059;                       % Vena cava elastance (mmHg/mL)
            R_sys = 1.0889 * 0.95;                  % Syst art resistance (mmHg*s/mL)
            % Heart valve paramenters
            R_mt = 0.0158;                          % Mitral valve resist (mmHg*s/mL)
            R_av = 0.018;                           % Aortic valve resist (mmHg*s/mL)
            R_tc = 0.0237;                          % Tricspd vlv resist (mmHg*s/mL)
            R_pv = 0.0055 * 1.5;                    % Pulmon vlv resist (mmHg*s/mL)
            % Heart failure param
            SVFact = 1.00;                          % Stress blood vol factor (uls)
            
            % Save all parameters in a sturcture to pass
            CVParam_Values = {period A B C E_es_lvf P_0_lvf ...
                lambda_lvf E_es_rvf P_0_rvf lambda_rvf E_es_pa ...
                E_es_pu R_pul E_es_ao E_es_vc R_sys R_mt R_av ...
                R_tc R_pv SVFact};
            CVParam_Fields = {'period' 'A' 'B' 'C' 'E_es_lvf' 'P_0_lvf' ...
                'lambda_lvf' 'E_es_rvf' 'P_0_rvf' 'lambda_rvf' ...
                'E_es_pa' 'E_es_pu' 'R_pul' 'E_es_ao' 'E_es_vc' ...
                'R_sys' 'R_mt' 'R_av' 'R_tc' 'R_pv' 'SVFact'};
            CVParam_Struct = cell2struct(CVParam_Values,CVParam_Fields,2);
        
        else
        
            % LOAD THE OPTIM PARAMS FROM AN EARLIER RUN
            % Left ventricle free wall parameters
            E_es_lvf = 2.8798;                      % LV free wall elast (mmHg/mL) 
            P_0_lvf = 0.0835; % 0.1203;             % LV ED pressure param (mmHg)
            lambda_lvf = 0.0295; % 0.033;           % LV ED pressure param (1/mL)
            % Right ventricle free wall parameters
            E_es_rvf = 0.585;                       % RV free wall elast (mmHg/mL) 
            P_0_rvf = 0.2157;                       % RV ED pressure param (mmHg)
            lambda_rvf = 0.023;                     % RV ED pressure param (1/mL)
            % Pulmonary artery and vein parameters
            E_es_pa = 0.369;                        % Pulm artery elastance (mmHg/mL)
            E_es_pu = 0.0073;                       % Pulm vein elastance (mmHg/mL)
            R_pul = 0.1552;                         % Pulm vasc resist (mmHg*s/mL)
            % Aortic and vena cava parameters
            E_es_ao = 1.1130; % 0.6913;             % Aorta elastance (mmHg/mL)
            E_es_vc = 0.0059;                       % Vena cava elastance (mmHg/mL)
            R_sys = 0.9687; % 1.0889;               % Syst art resistance (mmHg*s/mL)
            % Heart valve paramenters
            R_mt = 0.0158;                          % Mitral valve resist (mmHg*s/mL)
            R_av = 0.018;                           % Aortic valve resist (mmHg*s/mL)
            R_tc = 0.0237;                          % Tricspd vlv resist (mmHg*s/mL)
            R_pv = 0.0055;                          % Pulmon vlv resist (mmHg*s/mL)
            % Heart failure param
            SVFact = 2.4947; %1.00;                 % Stress blood vol factor (uls)
            
            % Save all parameters in a structure to pass
            CVParam_Values = {period A B C E_es_lvf P_0_lvf ...
                lambda_lvf E_es_rvf P_0_rvf lambda_rvf E_es_pa ...
                E_es_pu R_pul E_es_ao E_es_vc R_sys R_mt R_av ...
                R_tc R_pv SVFact};
            CVParam_Fields = {'period' 'A' 'B' 'C' 'E_es_lvf' 'P_0_lvf' ...
                'lambda_lvf' 'E_es_rvf' 'P_0_rvf' 'lambda_rvf' ...
                'E_es_pa' 'E_es_pu' 'R_pul' 'E_es_ao' 'E_es_vc' ...
                'R_sys' 'R_mt' 'R_av' 'R_tc' 'R_pv' 'SVFact'};
            CVParam_Struct = cell2struct(CVParam_Values,CVParam_Fields,2);
            
        end
        
    end

        
%% **********************************************************************************
%  Opt/Sim Params of    S M I T H   C A R D I O V A S C   M O D E L   S I M / O P T
% ***********************************************************************************
        
    % Pick the parameters we want to adjust with numbers shown below:
    
    % 1 --> E_es_lvf        LV free wall elast (mmHg/mL)        LEFT VENTRICLE
    % 2 --> P_0_lvf         LV ED pressure param (mmHg)
    % 3 --> lambda_lvf      LV ED pressure param (1/mL)
    % 4 --> E_es_rvf        RV free wall elast (mmHg/mL)        RIGHT VENTRICLE 
    % 5 --> P_0_rvf         RV ED pressure param (mmHg)
    % 6 --> lambda_rvf      RV ED pressure param (1/mL)
    % 7 --> E_es_pa         Pulm artery elastance (mmHg/mL)     PULMONARY 
    % 8 --> E_es_pu         Pulm vein elastance (mmHg/mL)       VASCULATURE
    % 9 --> R_pul           Pulm vasc resist (mmHg*s/mL)
    % 10 -> E_es_ao         Aorta elastance (mmHg/mL)           SYSTEMIC 
    % 11 -> E_es_vc         Vena cava elastance (mmHg/mL)       VASCULATURE
    % 12 -> R_sys           Syst art resistance (mmHg*s/mL)
    % 13 -> R_mt            Mitral valve resist (mmHg*s/mL)     HEART VALVES
    % 14 -> R_av            Aortic valve resist (mmHg*s/mL)
    % 15 -> R_tc            Tricspd vlv resist (mmHg*s/mL)
    % 16 -> R_pv            Pulmon vlv resist (mmHg*s/mL)
    % 17 -> SVFact          Stressed blood vol factor (uls)     CIRC BLOOD VOLUME
    
    % SELECT ADJUSTABLE PARAMETERS HERE
    AdjParams = [3;6;7;8;9;20;24];
    Num_AdjParams = size(AdjParams,1);
    
    % Bounds on all possible parameters
    LowBp_All(1) = log(0.100);      % LEFT          % LV free wall elast, E_es_lvf
    UpBp_All(1) = log(10.000);      %  VENTRICLE
    LowBp_All(2) = log(0.01);                       % LV ED pressure param, P_0_lvf
    UpBp_All(2) = log(5.00);
    LowBp_All(3) = log(0.005);                      % LV ED press param, lambda_lvf
    UpBp_All(3) = log(0.100);
    LowBp_All(4) = log(0.050);     % RIGHT         % RV free wall elast, E_es_rvf
    UpBp_All(4) = log(5.000);      %  VENTRICLE
    LowBp_All(5) = log(0.01);                      % RV ED pressure param, P_0_rvf
    UpBp_All(5) = log(5.00);
    LowBp_All(6) = log(0.005);                     % RV ED press param, lambda_rvf
    UpBp_All(6) = log(0.100);
    LowBp_All(7) = log(0.050);     % PULMONARY     % Pulm artery elastance, E_es_pa
    UpBp_All(7) = log(5.000);      %  VASCULATURE
    LowBp_All(8) = log(0.0005);                    % Pulm vein elastance, E_es_pu
    UpBp_All(8) = log(0.1000);
    LowBp_All(9) = log(0.005);                     % Pulm vasc resist, R_pul
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
            p(i) = log(CVParam_Values{ParamNum+4});
            AdjParam_Strngs(i,1) = CVParam_Fields(ParamNum+4);
            LowBp(i) = LowBp_All(ParamNum);
            UpBp(i) = UpBp_All(ParamNum);
        end
        
    else
        
        AdjParam_Strngs = {};
        
    end

    % Simulation time scales
    NumBeats_SS = 25;
    NumBeats_Sim = 5; 
    TSpan_SS = [0 NumBeats_SS * period];
    TSpan_Sim = [0 NumBeats_Sim * period];
    
    % Building the simulation parameter structure
    SimOptParam_Values = {AdjParam_Strngs TSpan_SS TSpan_Sim};
    SimOptParam_Fields = {'AdjParam_Strngs' 'TSpan_SS' 'TSpan_Sim'};
    SimOptParam_Struct = cell2struct(SimOptParam_Values, ...
        SimOptParam_Fields,2);
        
        
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
            fmc_Opts.MaxFunEvals = 5000;                % Set max num of fun evals
            fmc_Opts.MaxIter = 100;                     % Set max num of iterations
            if (Parallel_Flag == 1)
                fmc_Opts.UseParallel = 1;
            end
            
        % PERFORMING OPTIMIZATION
        if (Optim_Flag == 1)
            
            % fmincon optimization call
            p_Optim = fmincon(@Smith_ObjFun,p,[],[],[],[], ...
                LowBp,UpBp,[],fmc_Opts,PatData_Struct, ...
                CVParam_Struct,SimOptParam_Struct);
            
        elseif (Optim_Flag == 2 && Parallel_Flag == 0)
            
            % Genetic Algorithm Opimization - Serial computation
            ga_Opts = optimoptions('ga', ...
                'PopulationSize',250, ...
                'Display','iter', ...
                'MaxStallGenerations',10);
            Num_AdjParams = size(p,2);
            % Set objective function handle and 
            %  identify adjustable parameter vector, p
            Smith_ObjFun_Hndl = @(p) Smith_ObjFun(p, ...
                PatData_Struct,CVParam_Struct,SimOptParam_Struct);
            % Optimization call
            p_OptimGA = ga(Smith_ObjFun_Hndl,Num_AdjParams,[],[],[],[], ...
                LowBp,UpBp,[],ga_Opts);
            p_Optim = fmincon(@Smith_ObjFun,p_OptimGA,[],[],[],[], ...
                LowBp,UpBp,[],fmc_Opts,PatData_Struct, ...
                CVParam_Struct,SimOptParam_Struct);
            
        else
            
            % Genetic Algorithm Opimization - Parallel computation
            if (max(size(gcp)) == 0)                    % Parallel pool needed
                parpool                                 % Create the parallel pool
            end
            spmd
                warning('off','all');
            end  
            ga_Opts = optimoptions('ga', ...
                'PopulationSize',250, ...
                'Display','iter', ...
                'MaxStallGenerations',10, ...
                'UseParallel',true);
            Num_AdjParams = size(p,2);
            % Set objective function handle and 
            %  identify adjustable parameter vector, p
            Smith_ObjFun_Hndl = @(p) Smith_ObjFun(p, ...
                PatData_Struct,CVParam_Struct,SimOptParam_Struct);
            % Optimization call
            p_OptimGA = ga(Smith_ObjFun_Hndl,Num_AdjParams,[],[],[],[], ...
                LowBp,UpBp,[],ga_Opts);
            p_Optim = fmincon(@Smith_ObjFun,p_OptimGA,[],[],[],[], ...
                LowBp,UpBp,[],fmc_Opts,PatData_Struct, ...
                CVParam_Struct,SimOptParam_Struct);
            
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
%  Run Simulation of    S M I T H   C A R D I O V A S C   M O D E L   S I M / O P T
% ***********************************************************************************
           
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

    if (SmithParam_Flag == 1)
        % Calculating the initial condition on the volume due to the septum 
        %  deflection which is an implicit function requiring the use of fsolve
        time = 0;
        FSolve_Opts = optimset('Diagnostics','off', 'Display','off');
        SeptZF_Hndl = @(V_spt) SeptZF(V_spt,V_lv0,V_rv0,time,CVParam_Struct);
        [V_spt0,SeptZF_Val,ExitFlag,OutData] = fsolve(SeptZF_Hndl,-15, FSolve_Opts);
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
        X0(11) = V_spt0;
        
        % Build mass matrix for DAE solver
        M = eye(11);                                    % Put identity on diagonal
        M(11,11) = 0;                                   % Set last express as a ZFun
        % Set ODE/DAE options and time span
        ODE_Opts = odeset('Mass',M); 
        TSpan = [0 4];
        % Solve over the time span with ode15s
        [T_OutSS,X_OutSS] = ode15s(@dXdT_Smith, ...
            TSpan_SS,X0,ODE_Opts,CVParam_Struct);
        % Now solve over the simulation time span
        [T_Out,X_Out] = ode15s(@dXdT_Smith, ...
            TSpan_Sim,X_OutSS(end,:),ODE_Opts, ...
            CVParam_Struct);
        
    else
        
        X0(1) = V_lv0;
        X0(2) = V_rv0;
        X0(3) = V_pa0;
        X0(4) = V_pu0;
        X0(5) = V_ao0;
        X0(6) = V_vc0;
        
        % Solve over the time span with ode15s
        [T_OutSS,X_OutSS] = ode15s(@dXdT_SmithRed4, ...
            TSpan_SS,X0,[],CVParam_Struct);
        % Now solve over the simulation time span
        [T_Out,X_Out] = ode15s(@dXdT_SmithRed4, ...
            TSpan_Sim,X_OutSS(end,:),[],CVParam_Struct);
        
    end
    
    
%% **********************************************************************************
%  Plot Sim or Opt of   S M I T H   C A R D I O V A S C   M O D E L   S I M / O P T
% ***********************************************************************************

    % PLOT SOLUTION AND DATA
    %  In the case that we are not hand tuning we want to plot out the
    %  full final figure with data. If we are just plotting out a solution 
    %  Smith model results we can plot either with or without data
        
    % CALCULATING INTERMEDIATE PRESSURES TO PLOT
    Num_TOut = size(T_Out,1);                       % Number of time points
    P_LVSim = zeros(Num_TOut,1);                       % Preallocating matrices
    P_RVSim = zeros(Num_TOut,1);
    P_AOSim = zeros(Num_TOut,1);
    P_VCSim = zeros(Num_TOut,1);
    P_PASim = zeros(Num_TOut,1);
    P_PUSim = zeros(Num_TOut,1);
    V_LVESSim = 1000;
    V_LVEDSim = 0;
    % RUNNING MODEL TO GET INTERMEDIATES
    if (SmithParam_Flag == 1)
        
        for i = 1:Num_TOut
            VarOut = dXdT_Smith(T_Out(i),X_Out(i,:),CVParam_Struct,1);
            P_LVSim(i) = VarOut(1);
            P_RVSim(i) = VarOut(2);
            P_AOSim(i) = VarOut(3);
            P_VCSim(i) = VarOut(4);
            P_PASim(i) = VarOut(5);
            P_PUSim(i) = VarOut(6);
            V_LVESSim = min(V_LVESSim,X_Out(i,1));
            V_LVEDSim = max(V_LVEDSim,X_Out(i,1));
        end
        
    else
        
        for i = 1:Num_TOut
            VarOut = dXdT_SmithRed4(T_Out(i),X_Out(i,:),CVParam_Struct,1);
            P_LVSim(i) = VarOut(1);
            P_RVSim(i) = VarOut(2);
            P_AOSim(i) = VarOut(3);
            P_VCSim(i) = VarOut(4);
            P_PASim(i) = VarOut(5);
            P_PUSim(i) = VarOut(6);
            V_LVESSim = min(V_LVESSim,X_Out(i,1));
            V_LVEDSim = max(V_LVEDSim,X_Out(i,1));
        end
        
    end
        
    CO_Sim = ((V_LVEDSim - V_LVESSim) * Ave_HR) / 1000;
    if (HandTune_Flag == 0 && OptimBest_Flag == 0 && SmithParam_Flag == 0)
        Res_Optim = Smith_ObjFun(p_Optim,PatData_Struct, ...
            CVParam_Struct,SimOptParam_Struct);
    end
    
    % NOW GENERATE THE FIGURE 
    
    % Load previous Smith param model run if we are comparing runs
    if (Comp2Smith_Flag == 1)
        load SmithVars.mat
    end
   
    % Plot figure to check how close simulation 
    %  is to data or just simulation alone
    ScrSize = get(0,'ScreenSize');              % Getting screen size
    
    SixP_Fig = figure('Position', ...           % Positioning the figure
        [ScrSize(3)/50 ScrSize(4)/50 ...        %  on the screen
        ScrSize(3)/1.05 ScrSize(4)/1.05]); 
    
    if (SmithParam_Flag == 1)
        SupT_Hndl = suptitle({'Smith Params Normal'});
    elseif (Comp2Smith_Flag == 1)
        SupT_Hndl = suptitle({'Comparison to Smith Params'});
    else
        SupT_Hndl = suptitle({['D ID Number ' DID_Num]});
    end
    set(SupT_Hndl,'FontSize',24,'FontWeight','bold')
    
    % Subplot that compares measured and simulated pressures if data used
    subplot(2,3,1)
        plot(T_Out,P_RVSim,'-g', ...            % P_RV simulation
            'LineWidth',3, ...
            'DisplayName','P_{RV}')
        hold on
        if (Comp2Smith_Flag == 1)
            SmithP1 = plot(T_SimSmith, ...
                P_RVSimSmith,'--g', ...
                'LineWidth',3, ...
                'HandleVisibility','off');
            SmithP1.Color(4) = 0.25;
        end
        if (SmithParam_Flag == 0)
            plot([T_Out(1),T_Out(end)], ...     % P_RV systole data
                [P_RVsyst,P_RVsyst],'-.g', ...
                'LineWidth',1.5, ...
                'HandleVisibility','off')
            plot([T_Out(1),T_Out(end)], ...     % P_RV diastole data
                [P_RVdiast,P_RVdiast],':g', ...
                'LineWidth',1.5, ...
                'HandleVisibility','off')
        end
        plot(T_Out,P_AOSim,'-r', ...            % P_AO simulation
            'LineWidth',3, ...
            'DisplayName','P_{AO}')
        if (Comp2Smith_Flag == 1)
            SmithP2 = plot(T_SimSmith, ...
                P_AOSimSmith,'--r', ...
                'LineWidth',3, ...
                'HandleVisibility','off');
            SmithP2.Color(4) = 0.25;
        end
        if (SmithParam_Flag == 0)
            plot([T_Out(1),T_Out(end)], ...     % P_AO systole data
                [P_AOsyst,P_AOsyst],'-.r', ...
                'LineWidth',1.5, ...
                'HandleVisibility','off')
            plot([T_Out(1),T_Out(end)], ...     % P_AO diastole data
                [P_AOdiast,P_AOdiast],':r', ...
                'LineWidth',1.5, ...
                'HandleVisibility','off')
        end
        plot(T_Out,P_PASim,'-b', ...            % P_PA simulation
            'LineWidth',3, ...
            'DisplayName','P_{PA}')
        if (Comp2Smith_Flag == 1)
            SmithP3 = plot(T_SimSmith, ...
                P_PASimSmith,'--b', ...
                'LineWidth',3, ...
                'HandleVisibility','off');
            SmithP3.Color(4) = 0.25;
        end
        if (SmithParam_Flag == 0)
            plot([T_Out(1),T_Out(end)], ...     % P_PA systole data
                [P_PAsyst,P_PAsyst],'-.b', ...
                'LineWidth',1.5, ...
                'HandleVisibility','off')
            plot([T_Out(1),T_Out(end)], ...     % P_PA diastole data
                [P_PAdiast,P_PAdiast],':b', ...
                'LineWidth',1.5, ...
                'HandleVisibility','off')
        end
        plot(T_Out,P_PUSim,'-c', ...            % P_PU simulation
            'LineWidth',3, ...
            'DisplayName','P_{PU}')
        if (Comp2Smith_Flag == 1)
            SmithP4 = plot(T_SimSmith, ...
                P_PUSimSmith,'--c', ...
                'LineWidth',3, ...
                'HandleVisibility','off');
            SmithP4.Color(4) = 0.25;
        end
        if (SmithParam_Flag == 0)
            plot([T_Out(1),T_Out(end)], ...     % P_PCW average data
                [P_PCWave,P_PCWave],'-.c', ...
                'LineWidth',1.5, ...
                'HandleVisibility','off')
        end
        if (SmithParam_Flag == 0)
            PMax_SP1 = 1.45 * ...
                max(P_AOsyst,max(P_AOSim));
            text(0.5,0.95*PMax_SP1, ...         % CO data
                ['CO Data = ' num2str(CO_Thermo)])
            text(0.5,0.90*PMax_SP1, ...
                ['CO Sim  = ' num2str(CO_Sim)]) % CO simulation
        elseif (Comp2Smith_Flag == 1)
            PMax_SP1 = 1.45 * ...
                max([max(P_AOSim) max(P_AOSimSmith)]);
            text(0.5,0.95*PMax_SP1, ...         % CO Smith params
                ['CO Smith = ' num2str(CO_SimSmith)])
            text(0.5,0.90*PMax_SP1, ...
                ['CO Sim  = ' num2str(CO_Sim)]) % CO simulation
        else
            PMax_SP1 = 1.45 * max(P_AOSim);
            text(0.5,0.90*PMax_SP1, ...
                ['CO Sim  = ' num2str(CO_Sim)]) % CO simulation
        end
        xlim([0 5])                             % Formatting subplot
        ylim([-20 PMax_SP1])
        LegHndl1 = legend('show');
        set(LegHndl1,'Box','off','FontSize',8)
        set(gca,'FontSize',14,'FontWeight', ...
            'bold','Box','off')
        xlabel('Time (sec)','FontSize',20, ...
            'FontWeight','bold')
        ylabel('Pressures (mmHg)', ...
            'FontSize',20,'FontWeight','bold')
    
    % Subplot that shows simulated left and right ventricular volumes    
    subplot(2,3,2)
        plot(T_Out,X_Out(:,1),'-k', ...         % V_LV simulation    
            'LineWidth',3)
        hold on
        if (Comp2Smith_Flag == 1)
            SmithP5 = plot(T_SimSmith, ...
                V_LVSimSmith,'--k', ...
                'LineWidth',3, ...
                'HandleVisibility','off');
            SmithP5.Color(4) = 0.25;
        end
        plot(T_Out,X_Out(:,2),'-g', ...         % V_RV simulation 
            'LineWidth',3)
        if (Comp2Smith_Flag == 1)
            SmithP6 = plot(T_SimSmith, ...
                V_RVSimSmith,'--g', ...
                'LineWidth',3, ...
                'HandleVisibility','off');
            SmithP6.Color(4) = 0.25;
        end
        if (Comp2Smith_Flag == 1)
            VMax_SP2 = 1.20*max([max(X_Out(:,1)) ...
                max(X_Out(:,2)) ...
                max(V_LVSimSmith(:,1))]);
            VMin_SP2 = 0.85*min([min(X_Out(:,1)) ...
                min(X_Out(:,2)) ...
                min(V_LVSimSmith(:,1))]);
        else
            VMax_SP2 = 1.20*max(max(X_Out(:,1)), ...
                max(X_Out(:,2)));
            VMin_SP2 = 0.85*min(min(X_Out(:,1)), ...
                min(X_Out(:,2)));
        end
        xlim([0 5])                             % Formatting subplot
        ylim([VMin_SP2 VMax_SP2])
        LegHndl2 = legend('V_{LV}','V_{RV}');  
        set(LegHndl2,'Box','off','FontSize',8)
        set(gca,'FontSize',14, ...
            'FontWeight','bold','Box','off')
        xlabel('Time (sec)','FontSize',20, ...
            'FontWeight','bold')
        ylabel('Ventricular Volume (mL)', ...
            'FontSize',20,'FontWeight','bold')
        
    % Subplot that shows simulated left ventricular,
    %  aortic and pulmonary vein pressures
    subplot(2,3,4)
        plot(T_Out,P_LVSim,'-k','LineWidth',3)  % P_LV simulation
        hold on
        if (Comp2Smith_Flag == 1)
            SmithP7 = plot(T_SimSmith, ...
                P_LVSimSmith,'--k', ...
                'LineWidth',3, ...
                'HandleVisibility','off');
            SmithP7.Color(4) = 0.25;
        end
        plot(T_Out,P_AOSim,'-r','LineWidth',3)  % P_AO simulation
        if (Comp2Smith_Flag == 1)
            SmithP8 = plot(T_SimSmith, ...
                P_AOSimSmith,'--r', ...
                'LineWidth',3, ...
                'HandleVisibility','off');
            SmithP8.Color(4) = 0.25;
        end
        plot(T_Out,P_PUSim,'-c','LineWidth',3)  % P_PU simulation
        if (Comp2Smith_Flag == 1)
            SmithP9 = plot(T_SimSmith, ...
                P_PUSimSmith,'--c', ...
                'LineWidth',3, ...
                'HandleVisibility','off');
            SmithP9.Color(4) = 0.25;
        end
        PMax_SP4 = 1.30*max(max(P_LVSim), ...
            max(P_AOSim));
        xlim([0 5])                             % Formatting subplot
        ylim([-10 PMax_SP4])
        LegHndl4 = legend('P_{LV}', ...
            'P_{AO}','P_{PU}');
        set(LegHndl4,'Box','off','FontSize',8)
        set(gca,'FontSize',14, ...
            'FontWeight','bold','Box','off')
        xlabel('Time (sec)','FontSize',20, ...
            'FontWeight','bold')
        ylabel('Pressure (mmHg)', ...
            'FontSize',20,'FontWeight','bold')
    
    % Subplot that shows simulated right ventricular,
    %  pulmonary artery and vena cava pressures    
    subplot(2,3,5)
        plot(T_Out,P_RVSim,'-g','LineWidth',3)  % P_RV simulation
        hold on
        if (Comp2Smith_Flag == 1)
            SmithP10 = plot(T_SimSmith, ...
                P_RVSimSmith,'--g', ...
                'LineWidth',3, ...
                'HandleVisibility','off');
            SmithP10.Color(4) = 0.25;
        end
        plot(T_Out,P_PASim,'-b','LineWidth',3)  % P_PA simulation
        if (Comp2Smith_Flag == 1)
            SmithP11 = plot(T_SimSmith, ...
                P_PASimSmith,'--b', ...
                'LineWidth',3, ...
                'HandleVisibility','off');
            SmithP11.Color(4) = 0.25;
        end
        plot(T_Out,P_VCSim,'-m','LineWidth',3)  % P_VC simulation
        if (Comp2Smith_Flag == 1)
            SmithP12 = plot(T_SimSmith, ...
                P_VCSimSmith,'--m', ...
                'LineWidth',3, ...
                'HandleVisibility','off');
            SmithP12.Color(4) = 0.25;
        end
        PMax_SP5 = 1.30*max(max(P_RVSim), ...
            max(P_PASim));
        xlim([0 5])                             % Formatting subplot
        ylim([-5 PMax_SP5])
        LegHndl5 = legend('P_{RV}', ...
            'P_{PA}','P_{VC}');
        set(LegHndl5,'Box','off','FontSize',8)
        set(gca,'FontSize',14, ...
            'FontWeight','bold','Box','off')
        xlabel('Time (sec)','FontSize',20, ...
            'FontWeight','bold')
        ylabel('Pressure (mmHg)', ...
            'FontSize',20,'FontWeight','bold')
        
    % Subplot that shows simulated left and right ventricular pressure volume loops
    subplot(2,3,[3 6]) 
        LVVol_Max = max(X_Out(:,1));
        RVVol_Max = max(X_Out(:,2));
        if (SmithParam_Flag == 0)
            if (Comp2Smith_Flag == 1)
                LVVol_SmithMax = max(V_LVSimSmith);
                VMax_SP36 = max([LVVol_Max RVVol_Max ...
                    V_LVdiast LVVol_SmithMax]);
            else
                VMax_SP36 = max([LVVol_Max RVVol_Max V_LVdiast]);
            end
        else
            if (Comp2Smith_Flag == 1)
                LVVol_SmithMax = max(V_LVSimSmith);
                VMax_SP36 = max([LVVol_Max RVVol_Max LVVol_SmithMax]);
            else
                VMax_SP36 = max([LVVol_Max RVVol_Max]);
            end
        end
        if (Comp2Smith_Flag == 1)
            P_LVSmithMax = max(P_LVSimSmith);
            PMax_SP36 = max([max(P_LVSim) max(P_RVSim) P_LVSmithMax]);
        else
            PMax_SP36 = max(max(P_LVSim),max(P_RVSim));
        end
        PMin_SP36 = -10;
        plot(X_Out(:,1),P_LVSim,'-k', ...       % V_LV vs P_LV simulation
            'LineWidth',3)
        hold on
        if (Comp2Smith_Flag == 1)
            SmithP13 = plot(V_LVSimSmith, ...
                P_LVSimSmith,'--k', ...
                'LineWidth',3, ...
                'HandleVisibility','off');
            SmithP13.Color(4) = 0.15;
        end
        plot(X_Out(:,2),P_RVSim,'-g', ...       % V_RV vs P_RV simulation
            'LineWidth',3)
        if (Comp2Smith_Flag == 1)
            SmithP14 = plot(V_RVSimSmith, ...
                P_RVSimSmith,'--g', ...
                'LineWidth',3, ...
                'HandleVisibility','off');
            SmithP14.Color(4) = 0.15;
        end
        if ((SmithParam_Flag == 0) && (RHCEcho_Flag == 1))
            plot([V_LVsyst V_LVsyst], ...
                [0.5*PMin_SP36 1.05*PMax_SP36], ...
                '-.k','LineWidth',1.5)
            plot([V_LVdiast V_LVdiast], ...
                [0.5*PMin_SP36 1.05*PMax_SP36], ...
                ':k','LineWidth',1.5)
        end
        xlim([0 VMax_SP36*1.20])                % Formatting subplot
        ylim([PMin_SP36 1.15*PMax_SP36])
        LegHndl36 = legend('LV','RV');
        set(LegHndl36,'Box','off','FontSize',8)
        set(gca,'FontSize',14, ...
            'FontWeight','bold','Box','off')
        xlabel('LV Volume (mL)', ...
            'FontSize',20,'FontWeight','bold')
        ylabel('LV Pressure (mmHg)', ...
            'FontSize',20,'FontWeight','bold')
        
    if (SaveSmith_Flag == 1)
        
        T_SimSmith = T_Out;
        P_RVSimSmith = P_RVSim;
        P_AOSimSmith = P_AOSim;
        P_PASimSmith = P_PASim;
        P_PUSimSmith = P_PUSim;
        CO_SimSmith = CO_Sim;
        P_LVSimSmith = P_LVSim;
        P_VCSimSmith = P_VCSim;
        V_LVSimSmith = X_Out(:,1);
        V_RVSimSmith = X_Out(:,2);
        
        SaveFile1 = 'SmithVars.mat';
        save(SaveFile1, 'T_SimSmith', 'P_RVSimSmith', ...
            'P_AOSimSmith', 'P_PASimSmith', 'P_PUSimSmith', ...
            'CO_SimSmith', 'P_LVSimSmith', 'P_VCSimSmith',...
            'V_LVSimSmith', 'V_RVSimSmith')
        
    end
    if (HandTune_Flag == 0 && SmithParam_Flag == 0 && OptimBest_Flag == 0)
        Expp_Optim = exp(p_Optim);
        SaveFile2 = 'OptimParams.mat';
        save(SaveFile2,'Expp_Optim','AdjParam_Strngs','Res_Optim')
    end
    
    toc
