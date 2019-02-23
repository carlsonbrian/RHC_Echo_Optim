% ***********************************************************************************
%                    d X d T   O D E   F U N C T I O N   for 
%              R E D U C E D   S M I T H   C A R D I O V A S C U L A R   
%                S Y S T E M S   M O D E L   V E R S I O N   O N E
% ***********************************************************************************
%
%   This function contains the algebraic and differential expressions that describe
%   a reduced version of the Smith et al. cardiovascular system model (Med Eng Phys 
%   26:131, 2004) where the ventricular-ventricular interaction, valve inertances 
%   and all dead space/zero pressure volumes (except in pericardium) have been 
%   removed.
%
%   Model originally created on     19 September 2018
%   Model last modfied on           19 September 2018
%
%   Reproduced by       Brian Carlson
%                       Physiological Systems Dynamics Laboratory
%                       Department of Molecular and Integrative Physiology
%                       University of Michigan
%
% ***********************************************************************************
%  START OF  	   d X d t  for   R E D U C E D   S M I T H   M O D E L   V E R 1
% ***********************************************************************************

    function [Var_Out] = dXdT_SmithRed3(time,X,CVParam_Struct,varargin)

    % UNPACK FIXED PARAMETERS
    % Elastance function driver parameters
    period = CVParam_Struct.period;                 % Period of heart beat (s)
    A = CVParam_Struct.A;                           % Elastance function param (uls)
    B = CVParam_Struct.B;                           % Elastance fctn param (1/s^2)
    C = CVParam_Struct.C;                           % Elastance fctn param (s)
    % Pericardium calculation parameters
    P_0_pcd = CVParam_Struct.P_0_pcd;               % Unstress pericard press (kPa)
    lambda_pcd = CVParam_Struct.lambda_pcd;         % Pericard press exp term (1/mL)
    V_0_pcd = CVParam_Struct.V_0_pcd;               % Pericard zero P volume (mL)
    P_th = CVParam_Struct.P_th;                     % Thoracic pressure (kPa)
    % Left ventricle free wall parameters
    E_es_lvf = CVParam_Struct.E_es_lvf;             % LV free wall elastance (kPa/mL)
    P_0_lvf = CVParam_Struct.P_0_lvf;               % LV ED pressure param (kPa)
    lambda_lvf = CVParam_Struct.lambda_lvf;         % LV ED pressure param (1/mL)
    % Right ventricle free wall parameters
    E_es_rvf = CVParam_Struct.E_es_rvf;             % RV free wall elastance (kPa/mL)
    P_0_rvf = CVParam_Struct.P_0_rvf;               % RV ED pressure param (kPa)
    lambda_rvf = CVParam_Struct.lambda_rvf;         % RV ED pressure param (1/mL)
    % Pulmonary artery and vein parameters
    E_es_pa = CVParam_Struct.E_es_pa;               % Pulm artery elastance (kPa/mL)
    E_es_pu = CVParam_Struct.E_es_pu;               % Pulm vein elastance (kPa/mL)
    R_pul = CVParam_Struct.R_pul;                   % Pulm vasc resistance (kPa*s/mL)
    % Aortic and vena cava parameters
    E_es_ao = CVParam_Struct.E_es_ao;               % Aorta elastance (kPa/mL)
    E_es_vc = CVParam_Struct.E_es_vc;               % Vena cava elastance (kPa/mL)
    R_sys = CVParam_Struct.R_sys;                   % Syst art resistance (kPa*s/mL)
    % Heart valve paramenters
    R_mt = CVParam_Struct.R_mt;                     % Mitral valve resist (kPa*s/mL)
    R_av = CVParam_Struct.R_av;                     % Aortic valve resist (kPa*s/mL)
    R_tc = CVParam_Struct.R_tc;                     % Tricuspid vlv resist (kPa*s/mL)
    R_pv = CVParam_Struct.R_pv;                     % Pulmon vlv resist (kPa*s/mL)

    % Unpack the X vector
    V_lv = X(1);
    V_rv = X(2);
    V_pa = X(3);
    V_pu = X(4);
    V_ao = X(5);
    V_vc = X(6);
    
    %   <component name="driver_function">
    tau = time - (floor(time/period) * period);
	e_t = A * exp((-1) * B * (tau-C)^2);

	%   <component name="pericardium">
	V_pcd = V_lv + V_rv;
	P_pcd = P_0_pcd * (exp(lambda_pcd*(V_pcd-V_0_pcd))-1);
	P_peri = P_pcd + P_th;
    
    %   <component name="lvf_calculator">
    V_lvf = V_lv;
	P_es_lvf = E_es_lvf * V_lvf;
	P_ed_lvf = P_0_lvf * (exp(lambda_lvf * V_lvf) - 1);

	%   <component name="left_ventricle">
	P_lvf = (e_t * P_es_lvf) + ((1-e_t) * P_ed_lvf);
	P_lv = P_lvf + P_peri;
    P_pu = (E_es_pu * V_pu) + P_th;
    Q_mt = (P_pu - P_lv) / R_mt;
    P_ao = (E_es_ao * V_ao) + P_th;
    Q_av = (P_lv - P_ao) / R_av;
    if ((Q_mt < 0) && (Q_av < 0)) 
        dVlvdt = 0; 
    elseif (Q_mt < 0) 
        dVlvdt = (-1) * Q_av;
    elseif (Q_av < 0) 
        dVlvdt = Q_mt; 
    else
        dVlvdt = Q_mt - Q_av;
    end
    
    %   <component name="rvf_calculator">
    V_rvf = V_rv;
	P_es_rvf = E_es_rvf * V_rvf;
	P_ed_rvf = P_0_rvf * (exp(lambda_rvf * V_rvf) - 1);
    
	%   <component name="right_ventricle">
	P_rvf = (e_t * P_es_rvf) + ((1-e_t) * P_ed_rvf);
	P_rv = P_rvf + P_peri;
    P_vc = (E_es_vc * V_vc) + P_th;
    Q_tc = (P_vc - P_rv) / R_tc;
    P_pa = (E_es_pa * V_pa) + P_th;
    Q_pv = (P_rv - P_pa) / R_pv;
    if ((Q_tc < 0) && (Q_pv < 0)) 
        dVrvdt = 0; 
    elseif (Q_tc < 0) 
        dVrvdt = (-1) * Q_pv;
    elseif (Q_pv < 0) 
        dVrvdt = Q_tc; 
    else
        dVrvdt = Q_tc - Q_pv;
    end

	%   <component name="pulmonary_artery and vein">
	Q_pul = (P_pa - P_pu) / R_pul;
    if (Q_pv < 0)
        dVpadt = -1 * Q_pul;
    else
        dVpadt = Q_pv - Q_pul;
    end
    if (Q_mt < 0)
        dVpudt = Q_pul;
    else
        dVpudt = Q_pul - Q_mt;
    end
    
	%   <component name="aorta and vena cava">
    Q_sys = (P_ao - P_vc) / R_sys;
    if (Q_av < 0)
        dVaodt = -1 * Q_sys;
    else
        dVaodt = Q_av - Q_sys;
    end
    if (Q_tc < 0)
        dVvcdt = Q_sys;
    else
        dVvcdt = Q_sys - Q_tc;
    end    
    
    if (isempty(varargin))       
        RoC(1) = dVlvdt;
        RoC(2) = dVrvdt;
        RoC(3) = dVpadt;
        RoC(4) = dVpudt;
        RoC(5) = dVaodt;
        RoC(6) = dVvcdt;
        Var_Out = RoC';     
    else     
        CalcVars(1) = P_lv;
        CalcVars(2) = P_rv;
        CalcVars(3) = P_ao;
        CalcVars(4) = P_vc;
        CalcVars(5) = P_pa;
        CalcVars(6) = P_pu;
        Var_Out = CalcVars';
    end
end

