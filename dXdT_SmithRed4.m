% ***********************************************************************************
%                    d X d T   O D E   F U N C T I O N   for 
%              R E D U C E D   S M I T H   C A R D I O V A S C U L A R   
%                S Y S T E M S   M O D E L   V E R S I O N   O N E
% ***********************************************************************************
%
%   This function contains the algebraic and differential expressions that describe
%   a reduced version of the Smith et al. cardiovascular system model (Med Eng Phys 
%   26:131, 2004) where the ventricular-ventricular interaction, valve inertances, 
%   all dead space/zero pressure volumes, pericardium and thoracic chambers have 
%   been removed.
%
%   Model originally created on     20 September 2018
%   Model last modfied on           24      July 2020
%
%   Reproduced by       Brian Carlson
%                       Physiological Systems Dynamics Laboratory
%                       Department of Molecular and Integrative Physiology
%                       University of Michigan
%
% ***********************************************************************************
%  START OF  	   d X d t  for   R E D U C E D   S M I T H   M O D E L   V E R 1
% ***********************************************************************************

    function [Var_Out] = dXdT_SmithRed4(time,X,DriverP_Struct, ...
        CVParam_Struct,varargin)

    % UNPACK FIXED PARAMETERS
    % Elastance function driver parameters
    period = DriverP_Struct.period;                 % Period of heart beat (s)
    A = CVParam_Struct.A;                           % Elastance function param (uls)
    B = DriverP_Struct.B;                           % Elastance fctn param (1/s^2)
    C = DriverP_Struct.C;                           % Elastance fctn param (s)
    % Left ventricle free wall parameters
    E_es_lvf = CVParam_Struct.E_es_lvf;             % LV free wall elstnce (mmHg/mL)
    P_0_lvf = CVParam_Struct.P_0_lvf;               % LV ED pressure param (mmHg)
    lambda_lvf = CVParam_Struct.lambda_lvf;         % LV ED pressure param (1/mL)
    % Right ventricle free wall parameters
    E_es_rvf = CVParam_Struct.E_es_rvf;             % RV free wall elstnce (mmHg/mL)
    P_0_rvf = CVParam_Struct.P_0_rvf;               % RV ED pressure param (mmHg)
    lambda_rvf = CVParam_Struct.lambda_rvf;         % RV ED pressure param (1/mL)
    % Pulmonary artery and vein parameters
    E_es_pa = CVParam_Struct.E_es_pa;               % Pulm arterial elstnce (mmHg/mL)
    E_es_pu = CVParam_Struct.E_es_pu;               % Pulm venous elastance (mmHg/mL)
    R_pul = CVParam_Struct.R_pul;                   % Pulm vasc rsistnce (mmHg*s/mL)
    % Aortic and vena cava parameters
    E_es_sa = CVParam_Struct.E_es_sa;               % Syst arter elastance (mmHg/mL)
    E_es_sv = CVParam_Struct.E_es_sv;               % Syst venous elastance (mmHg/mL)
    R_sys = CVParam_Struct.R_sys;                   % Syst vasc rsistnce (mmHg*s/mL)
    % Heart valve paramenters
    R_mt = CVParam_Struct.R_mt;                     % Mitral valve resist (mmHg*s/mL)
    R_av = CVParam_Struct.R_av;                     % Aortic valve resist (mmHg*s/mL)
    R_tc = CVParam_Struct.R_tc;                     % Tricspd vlv resist (mmHg*s/mL)
    R_pv = CVParam_Struct.R_pv;                     % Pulmon vlv resist (mmHg*s/mL)

    % Unpack the X vector
    V_lv = X(1);
    V_rv = X(2);
    V_pa = X(3);
    V_pu = X(4);
    V_sa = X(5);
    V_sv = X(6);
    
    % Elastance function driving heart systole and diastole
    tau = time - (floor(time/period) * period);
	e_t = A * exp((-1) * B * (tau-C)^2);
    
    % Left ventricular free wall end systole and end diastole PVRs
    V_lvf = V_lv;
	P_es_lvf = E_es_lvf * V_lvf;
	P_ed_lvf = P_0_lvf * (exp(lambda_lvf * V_lvf) - 1);

	% Left ventricular pressure development
	P_lvf = (e_t * P_es_lvf) + ((1-e_t) * P_ed_lvf);
	P_lv = P_lvf;
    
    % Pulmonary venous pressure development 
    %  and flow through mitral valve
    P_pu = E_es_pu * V_pu;
    Q_mt = (P_pu - P_lv) / R_mt;
    
    % Systemic arterial pressure development 
    %  and flow through aortic valve
    P_sa = E_es_sa * V_sa;
    Q_av = (P_lv - P_sa) / R_av;
    
    % Change in left ventricular volume based on 
    %  mitral and aortic valve flow
    if ((Q_mt <= 0) && (Q_av <= 0)) 
        dVlvdt = 0; 
    elseif (Q_mt <= 0) 
        dVlvdt = (-1) * Q_av;
    elseif (Q_av <= 0) 
        dVlvdt = Q_mt; 
    else
        dVlvdt = Q_mt - Q_av;
    end
    
    % Right ventricular free wall end systole and end diastole PVRs
    V_rvf = V_rv;
	P_es_rvf = E_es_rvf * V_rvf;
	P_ed_rvf = P_0_rvf * (exp(lambda_rvf * V_rvf) - 1);
    
	% Left ventricular pressure development
	P_rvf = (e_t * P_es_rvf) + ((1-e_t) * P_ed_rvf);
	P_rv = P_rvf;
    
    % Systemic venous pressure development 
    %  and flow through tricuspid valve
    P_sv = E_es_sv * V_sv;
    Q_tc = (P_sv - P_rv) / R_tc;
    
    % Pulmonary arterial pressure development 
    %  and flow through pulmonary valve
    P_pa = E_es_pa * V_pa;
    Q_pv = (P_rv - P_pa) / R_pv;
    
    % Change in right ventricular volume based on 
    %  tricuspid and pulmonary valve flow
    if ((Q_tc <= 0) && (Q_pv <= 0)) 
        dVrvdt = 0; 
    elseif (Q_tc <= 0) 
        dVrvdt = (-1) * Q_pv;
    elseif (Q_pv <= 0) 
        dVrvdt = Q_tc; 
    else
        dVrvdt = Q_tc - Q_pv;
    end

	% Change in pulmonary vasculature volume based on
    %  pulmonary and mitral valve flow
	Q_pul = (P_pa - P_pu) / R_pul;
    if (Q_pv <= 0)
        dVpadt = -1 * Q_pul;
    else
        dVpadt = Q_pv - Q_pul;
    end
    if (Q_mt <= 0)
        dVpudt = Q_pul;
    else
        dVpudt = Q_pul - Q_mt;
    end
    
	% Change in systemic vasculature volume based on
    %  aortic and tricuspid valve flow
    Q_sys = (P_sa - P_sv) / R_sys;
    if (Q_av <= 0)
        dVsadt = -1 * Q_sys;
    else
        dVsadt = Q_av - Q_sys;
    end
    if (Q_tc <= 0)
        dVsvdt = Q_sys;
    else
        dVsvdt = Q_sys - Q_tc;
    end    
    
    % Returns either the derivatives during integration or the
    %  intermediate pressures during post integration runs
    if (isempty(varargin))       
        RoC(1) = dVlvdt;
        RoC(2) = dVrvdt;
        RoC(3) = dVpadt;
        RoC(4) = dVpudt;
        RoC(5) = dVsadt;
        RoC(6) = dVsvdt;
        Var_Out = RoC';     
    else     
        CalcVars(1) = P_lv;
        CalcVars(2) = P_rv;
        CalcVars(3) = P_sa;
        CalcVars(4) = P_sv;
        CalcVars(5) = P_pa;
        CalcVars(6) = P_pu;
        Var_Out = CalcVars';
    end
end

