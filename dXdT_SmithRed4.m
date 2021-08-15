% ***********************************************************************************
%                    d X d T   O D E   F U N C T I O N   for 
%              R E D U C E D   S M I T H   C A R D I O V A S C U L A R   
%                S Y S T E M S   M O D E L   V E R S I O N   F O U R
% ***********************************************************************************
%
%   This function contains the algebraic and differential expressions that describe
%   a reduced version of the Smith et al. cardiovascular system model (Med Eng Phys 
%   26:131, 2004) where the ventricular-ventricular interaction, valve inertances, 
%   all dead space/zero pressure volumes, pericardium and thoracic chambers have 
%   been removed.
%
%   Model originally created on     20 September 2018
%   Model last modfied on           21      July 2021
%
%   Reproduced by       Brian Carlson
%                       Physiological Systems Dynamics Laboratory
%                       Department of Molecular and Integrative Physiology
%                       University of Michigan
%
% ***********************************************************************************
%  START OF  	   d X d t  for   R E D U C E D   S M I T H   M O D E L   V E R 4
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
    E_lv = CVParam_Struct.E_lv;                     % LV free wall elstnce (mmHg/mL)
    P_0lv = CVParam_Struct.P_0lv;                   % LV ED pressure param (mmHg)
    lambda_lv = CVParam_Struct.lambda_lv;           % LV ED pressure param (1/mL)
    % Right ventricle free wall parameters
    E_rv = CVParam_Struct.E_rv;                     % RV free wall elstnce (mmHg/mL)
    P_0rv = CVParam_Struct.P_0rv;                   % RV ED pressure param (mmHg)
    lambda_rv = CVParam_Struct.lambda_rv;           % RV ED pressure param (1/mL)
    % Pulmonary artery and vein parameters
    E_pa = CVParam_Struct.E_pa;                     % Pulm arterial elstnce (mmHg/mL)
    E_pv = CVParam_Struct.E_pv;                     % Pulm venous elastance (mmHg/mL)
    R_pul = CVParam_Struct.R_pul;                   % Pulm vasc rsistnce (mmHg*s/mL)
    % Aortic and vena cava parameters
    E_sa = CVParam_Struct.E_sa;                     % Syst arter elastance (mmHg/mL)
    E_sv = CVParam_Struct.E_sv;                     % Syst venous elastance (mmHg/mL)
    R_sys = CVParam_Struct.R_sys;                   % Syst vasc rsistnce (mmHg*s/mL)
    % Heart valve paramenters
    R_mval = CVParam_Struct.R_mval;                 % Mitral valve resist (mmHg*s/mL)
    R_aval = CVParam_Struct.R_aval;                 % Aortic valve resist (mmHg*s/mL)
    R_tval = CVParam_Struct.R_tval;                 % Tricspd vlv resist (mmHg*s/mL)
    R_pval = CVParam_Struct.R_pval;                 % Pulmon vlv resist (mmHg*s/mL)

    % Unpack the X vector
    V_lv = X(1);
    V_rv = X(2);
    V_pa = X(3);
    V_pv = X(4);
    V_sa = X(5);
    V_sv = X(6);
    
    % Elastance function driving heart systole and diastole
    tau = time - (floor(time/period) * period);
	e_t = A * exp((-1) * B * (tau-C)^2);
    
    % Left ventricular end systole and end diastole PVRs
	P_es_lv = E_lv * V_lv;
	P_ed_lv = P_0lv * (exp(lambda_lv * V_lv) - 1);

	% Left ventricular pressure development
	P_lv = (e_t * P_es_lv) + ((1-e_t) * P_ed_lv);

    % Pulmonary venous pressure development 
    %  and flow through mitral valve
    P_pv = E_pv * V_pv;
    Q_mval = (P_pv - P_lv) / R_mval;
    
    % Systemic arterial pressure development 
    %  and flow through aortic valve
    P_sa = E_sa * V_sa;
    Q_aval = (P_lv - P_sa) / R_aval;
    
    % Change in left ventricular volume based on 
    %  mitral and aortic valve flow
    if ((Q_mval <= 0) && (Q_aval <= 0)) 
        dVlvdt = 0; 
    elseif (Q_mval <= 0) 
        dVlvdt = (-1) * Q_aval;
    elseif (Q_aval <= 0) 
        dVlvdt = Q_mval; 
    else
        dVlvdt = Q_mval - Q_aval;
    end
    
    % Right ventricular free wall end systole and end diastole PVRs
	P_es_rv = E_rv * V_rv;
	P_ed_rv = P_0rv * (exp(lambda_rv * V_rv) - 1);
    
	% Right ventricular pressure development
	P_rv = (e_t * P_es_rv) + ((1-e_t) * P_ed_rv);
    
    % Systemic venous pressure development 
    %  and flow through tricuspid valve
    P_sv = E_sv * V_sv;
    Q_tval = (P_sv - P_rv) / R_tval;
    
    % Pulmonary arterial pressure development 
    %  and flow through pulmonary valve
    P_pa = E_pa * V_pa;
    Q_pval = (P_rv - P_pa) / R_pval;
    
    % Change in right ventricular volume based on 
    %  tricuspid and pulmonary valve flow
    if ((Q_tval <= 0) && (Q_pval <= 0)) 
        dVrvdt = 0; 
    elseif (Q_tval <= 0) 
        dVrvdt = (-1) * Q_pval;
    elseif (Q_pval <= 0) 
        dVrvdt = Q_tval; 
    else
        dVrvdt = Q_tval - Q_pval;
    end

	% Change in pulmonary vasculature volume based on
    %  pulmonary and mitral valve flow
	Q_pul = (P_pa - P_pv) / R_pul;
    if (Q_pval <= 0)
        dVpadt = -1 * Q_pul;
    else
        dVpadt = Q_pval - Q_pul;
    end
    if (Q_mval <= 0)
        dVpvdt = Q_pul;
    else
        dVpvdt = Q_pul - Q_mval;
    end
    
	% Change in systemic vasculature volume based on
    %  aortic and tricuspid valve flow
    Q_sys = (P_sa - P_sv) / R_sys;
    if (Q_aval <= 0)
        dVsadt = -1 * Q_sys;
    else
        dVsadt = Q_aval - Q_sys;
    end
    if (Q_tval <= 0)
        dVsvdt = Q_sys;
    else
        dVsvdt = Q_sys - Q_tval;
    end    
    
    % Returns either the derivatives during integration or the
    %  intermediate pressures during post integration runs
    if (isempty(varargin))       
        RoC(1) = dVlvdt;
        RoC(2) = dVrvdt;
        RoC(3) = dVpadt;
        RoC(4) = dVpvdt;
        RoC(5) = dVsadt;
        RoC(6) = dVsvdt;
        Var_Out = RoC';     
    else     
        CalcVars(1) = P_lv;
        CalcVars(2) = P_rv;
        CalcVars(3) = P_sa;
        CalcVars(4) = P_sv;
        CalcVars(5) = P_pa;
        CalcVars(6) = P_pv;
        Var_Out = CalcVars';
    end
end

