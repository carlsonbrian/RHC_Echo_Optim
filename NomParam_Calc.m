%% **********************************************************************************
%              N O M I N A L   P A R A M E T E R   C A L C U L A T O R   
% ***********************************************************************************
%
%   This function calculates the nominal parameter values used to set the fixed
%   parameters in a patient specific manner and provide the initial estimated
%   cardiovascular function for visualization or starting point of a gradient-based
%   optimization method. We first start by estimating the volume distribution 
%   across compartments in this model and then the percentages of stressed volume 
%   in each compartment such that the stressed blood volume is what is 
%   prescribed in the nominal parameter calculation. The nominal values of volume
%   are calculated in accordance with the Beneken chapter (Beneken, J.E.W., A 
%   physical approach to hemodynamic aspects of the human cardiovascular system. 
%   Physical bases of circulatory transport: regulation and exchange. W.B. 
%   Saunders, Philadelphia, PA, USA., 1967).  The difference here is that we can 
%   have the total stressed blood volume vary from the Beneken value of 18.75% and
%   the code will recalculate an initial stressed volume distribution that is 
%   appropriate. For example, in many of our codes we have taken from literature
%   that the stressed volume is actually 30% of the total blood volume so the 
%   initial volume in each compartment will be larger than given in the table below.
%
%   Here is a summary of Table 1-1 values we use as a starting point from Beneken:
%
%      Compartment      Str Vol (mL)    Unstr Vol (mL)     Total Vol (mL)
%      -----------     -------------    --------------     --------------
%   Systemic arteries       160              425                 585
%   Systemic veins          219             2697                2916
%   Pulmonary arteries       69               50                 119
%   Pulmonary veins          54              460                 514
%   Left ventricle          125                0                 125
%   Left atrium              50               30                  80
%   Right ventricle         125                0                 125
%   Right atrium             50               30                  80
%   ----------------------------------------------------------------------
%   Totals                  852             3692                4544
%
%   Once the nominal volumes have been calculated then we can use them in the
%   calculations of elastance and the remaining parameters can be nominally
%   calculated using the clinical measures. One parameter in the end diastolic
%   pressure-volume relationship for the left and the right ventricle cannot be 
%   calculated - the pressure scaling terms, P_0lvf and P_0rvf. These are set to
%   the values used in the original Smith et al. model.
%
%   Model originally created on     26  June 2020
%   Model last modfied on            8 March 2021
%
%   Developed by        Brian Carlson
%                       Physiological Systems Dynamics Laboratory
%                       Department of Molecular and Integrative Physiology
%                       University of Michigan
%
%% ***********************************************************************************
%  START OF        N O M I N A L   B L O O D   V O L U M E   C A L C U L A T O R    
% ***********************************************************************************

    function NomParam_Struct = NomParam_Calc(PatData_Struct, ...
        CVParam_Struct,RHCData_Struct,EchoData_Struct)
    
    
    %% GET PATIENT VALUES NEEDED
    % Unpack patient data
    Hgt = PatData_Struct.Hgt;
    BW = PatData_Struct.BW;
    Sex = PatData_Struct.Sex;
    % Unpack parameters fixed for the purpose of nominal 
    %  parameter calculations based on normal cardiovascular 
    %  values or set to reflect heart failure conditions
    P_0_lvf = CVParam_Struct.P_0_lvf;
    P_0_rvf = CVParam_Struct.P_0_rvf;
    SVFact = CVParam_Struct.SVFact;
    % Unpack patient measures
    P_RVsyst = RHCData_Struct.P_RVsyst;
    P_RVdiast = RHCData_Struct.P_RVdiast;
    P_PAsyst = RHCData_Struct.P_PAsyst;
    P_PAdiast = RHCData_Struct.P_PAdiast;
    P_PCWave = RHCData_Struct.P_PCWave;
    P_AOsyst = RHCData_Struct.P_AOsyst;
    P_AOdiast = RHCData_Struct.P_AOdiast;
    CO_Fick = RHCData_Struct.CO_Fick;
    CO_Thermo = RHCData_Struct.CO_Thermo;
    V_LVsyst = EchoData_Struct.V_LVsyst;
    V_LVdiast = EchoData_Struct.V_LVdiast;
    ID_LVsyst = EchoData_Struct.ID_LVsyst;  
    ID_LVdiast = EchoData_Struct.ID_LVdiast;  
    % Selecting the best available cardiac output for RHC
    if (CO_Thermo == -1)
        CO_RHC = CO_Fick;
    else
        CO_RHC = CO_Thermo;
    end
    
    % If the LV volume is not calculated with 2D Simpson's method
    %  then use Teichholz formula with the 1D top of the ventricle
    %  diameter to estimate LV volume in systole
    if (V_LVsyst == -1)
        % Calculate the LV systolic volume using the Teichholz expression
        V_LVsyst = ((ID_LVsyst/10)^3) / ...       
            ((6/pi) * ((0.075 * (ID_LVsyst/10)) + 0.18));
        V_LVdiast = ((ID_LVdiast/10)^3) / ...       
            ((6/pi) * ((0.075 * (ID_LVdiast/10)) + 0.18));
    end
        
    
    %% RECALCULATE STRESSED VOLUME FRACTIONS FROM BENEKEN FOR 0.30 * SVFACT
    % Values from Table 1-1 in Beneken chapter (some values not used but
    %  are inserted here for completeness)
    VsB_sa = 160;           VuB_sa = 425;           VtB_sa = 585;
    VsB_sv = 219;           VuB_sv = 2697;          VtB_sv = 2916;
    VsB_pa = 69;            VuB_pa = 50;            VtB_pa = 119;
    VsB_pu = 54;            VuB_pu = 460;           VtB_pu = 514;
    VsB_lv = 125;           VuB_lv = 0;             VtB_lv = 125;
    VsB_la = 50;            VuB_la = 30;            VtB_la = 80;
    VsB_rv = 125;           VuB_rv = 0;             VtB_rv = 125;
    VsB_ra = 50;            VuB_ra = 30;            VtB_ra = 80;
    VsB_tot = 852;          VuB_tot = 3692;         VtB_tot = 4544;
    
    % The stressed volume in the left and right heart in Beneken assumes
    %  a 100% ejection fraction in the ventricles and a 60% ejection fraction
    %  in the atria so assuming a more realistic 70% and 50% ejection fraction
    %  in the ventricles and atria respectively gives
    VsBNew_lv = 0.70 * VtB_lv;
    VsBNew_la = 0.50 * VtB_la;
    VsBNew_rv = 0.70 * VtB_rv;
    VsBNew_ra = 0.50 * VtB_ra;
    
    % Recruiting volume over the Beneken values will come from only the
    %  systemic and pulmonary circulations and not from the heart so we 
    %  first take the stressed volume directly from Beneken and then
    %  subtract off the heart chamber stressed volumes
    VsB_nh = VsB_tot - VsB_lv - VsB_la - VsB_rv - VsB_ra;
    
    % But now assume we have the stressed volume fraction specified
    %  in this patient with the now more realistic stressed volumes 
    %  calculated above but with the Beneken total blood volume
    VsBP_tot = 0.30 * SVFact * VtB_tot;
    VsBP_nh = VsBP_tot - VsBNew_lv - VsBNew_la - VsBNew_rv - VsBNew_ra;
    
    % The difference between these two stressed volumes is the amount
    %  that we have to recruit over and above the stressed volumes given
    %  in Beneken Table 1-1
    VsRec_tot = VsBP_nh - VsB_nh;
    
    % Now we have to decide where to recruit this from. If we take a 
    %  straight percentage based on total volume in each compartment
    %  we run the risk of recruiting all of the volume in the compartments
    %  that have little left to give (e.g. pulmonary arteries). So this
    %  should be calculated based on the fraction of the unstressed volume 
    %  in each compartment with respect to the total unstressed volume,
    %  VuB_tot. This way all compartments would approach 0% unstressed
    %  volume equally as the stressed volume appoached 100% 
    VsBRec_sa = VsRec_tot * (VuB_sa / VuB_tot);
    VsBRec_sv = VsRec_tot * (VuB_sv / VuB_tot);
    VsBRec_pa = VsRec_tot * (VuB_pa / VuB_tot);
    VsBRec_pu = VsRec_tot * (VuB_pu / VuB_tot);
    
    % So adding these recruited volumes to the Beneken Table 1-1 values
    %  assumed at 18.75% stressed blood volume will give the volumes in
    %  in the Beneken example with 0.30 * SVFact total stressed volume.
    VsBNew_sa = VsB_sa + VsBRec_sa;
    VsBNew_sv = VsB_sv + VsBRec_sv;
    VsBNew_pa = VsB_pa + VsBRec_pa;
    VsBNew_pu = VsB_pu + VsBRec_pu;
    
    % Dividing these new volumes by the total compartment volumes from 
    %  Benenken Table 1-1 gives the fraction of the total compartment
    %  volume that is stressed with a 0.30 * SVFact total stressed volume.
    VsBFrac_sa = VsBNew_sa / VtB_sa;
    VsBFrac_sv = VsBNew_sv / VtB_sv;
    VsBFrac_pa = VsBNew_pa / VtB_pa;
    VsBFrac_pu = VsBNew_pu / VtB_pu;
    

    %% NOW TRANSLATED THE RECALCULATED STRESSED VOLUME FRACTIONS TO PATIENT
    % To calculate the stressed blood volume first the  
    %  total blood volume must be estimated from height weight and sex. 
    if (Sex == 'M')
        TotBV = ((0.3669 * (Hgt/100)^3) + ...
            (0.03219 * BW) + 0.6041) * 1000;
    else
        TotBV = ((0.3561 * (Hgt/100)^3) + ...
            (0.03308 * BW) + 0.1833) * 1000;
    end
    StrBV = SVFact * 0.30 * TotBV;
    
    % Get the total blood volume in each compartment using the values from
    %  the Beneken chapter Table 1-1
    Vsa_tot = StrBV * (VtB_sa / VtB_tot);
    Vsv_tot = StrBV * (VtB_sv / VtB_tot);
    Vpa_tot = StrBV * (VtB_pa / VtB_tot);
    Vpu_tot = StrBV * (VtB_pu / VtB_tot);
    
    % Now calculate the stressed volumes using the 
    %  newly calculated fractions from the previous section
    Vs_sa = VsBFrac_sa * Vsa_tot;
    Vs_sv = VsBFrac_sv * Vsv_tot;
    Vs_pa = VsBFrac_pa * Vpa_tot;
    Vs_pu = VsBFrac_pu * Vpu_tot;
    
    % Returning these values using the variable names used in the main script
    V_SAsyst = Vs_sa;
    V_SVsyst = Vs_sv;
    V_PAsyst = Vs_pa;
    V_PUsyst = Vs_pu;
    
    
    %% CALCULATING THE NOMINAL PARAMETERS
    % LEFT VENTRICLE FREE WALL PARAMETERS
    % Calculate nominal E_es_lvf
    PLVsyst_Nom = 1.025 * P_AOsyst;
    Eeslvf_Nom = PLVsyst_Nom / V_LVsyst;                % Calculated from data
    % Set P_0_lvf to normal Smith value
    P0lvf_Nom = P_0_lvf;                                % Normal value
    % Calculate nominal lambda_lvf
    PPUpp_Nom = 0.20 * (P_PAsyst - P_PAdiast);          % Nom pulm venous pulse P
    PPUdiast_Nom = P_PCWave - (PPUpp_Nom / 3);          % Nom pulm venous diast P
    PLVdiast_Nom = 0.975 * PPUdiast_Nom;
    lambdalvf_Nom = (log(PLVdiast_Nom) - ...
        log(P0lvf_Nom)) / V_LVdiast;
    
    % RIGHT VENTRICULAR FREE WALL PARAMETERS
    % Calculate nominal E_es_rvf
    VRVsyst_Nom = 0.90 * V_LVsyst;                      % Used in Heart Tx paper
    Eesrvf_Nom = P_RVsyst / VRVsyst_Nom;                % Calculated from data
    % Set P_0_rvf to normal Smith value
    P0rvf_Nom = P_0_rvf;                                % Normal value
    % Calculate nominal lambda_rvf
    VRVdiast_Nom = 0.90 * V_LVdiast;                    % Used in Heart Tx paper
    PSVpp_Nom = 0.05 * (P_AOsyst - P_AOdiast);
    % Calculating the nominal value of lambda_rvf 
    if (P_RVdiast <= 0)                                 % This takes care of when
        PRVdiast_Nom = PSVpp_Nom/2;                     %  P_RVdiast is zero (or
        lambdarvf_Nom = (log(PRVdiast_Nom) - ...        %  negative) in the
            log(P0rvf_Nom)) / VRVdiast_Nom;             %  clinical data
    else
        lambdarvf_Nom = (log(P_RVdiast) - ...
            log(P0rvf_Nom)) / VRVdiast_Nom;  
    end
    
    % Pulmonary artery and vein parameters
    Eespa_Nom = (P_PAsyst - P_PAdiast) / V_PAsyst;      % Nom pulm art elstnce
    Eespu_Nom = PPUpp_Nom / V_PUsyst;                   % Nom pulm venous elstnce
    % The factor in the second nominal value of P_PUsyst is based on the 
    %  average P_PUsyst/P_PAsyst for all patients in HFpEF/HFrEF study where 
    %  P_PCWave is less than P_PAsyst. Patient BCHF14 has P_PCWave = 78 mmHg and 
    %  P_PAsyst = 65 mmHg so this catches this type of mismeasurement. In this case
    %  the P_PUsyst from P_PAsyst based on the average of the other 30 patients in
    %  this study. The factor (0.4854) representing this ratio of P_PUsyst/P_PAsyst 
    %  might be adjusted as our patient pool grows.
    P_PAave = (P_PAsyst/3) + ((2*P_PAdiast)/3);
    PPUsyst_Nom1 = P_PCWave + ((2*PPUpp_Nom)/3);        % Nom pulm venous systolic P
    PPUsyst_Nom2 = 0.4854 * P_PAsyst;                   % If P_PCWave is too big
    if ((0.975 * P_PAave) >= PPUsyst_Nom1)
        Rpul_Nom = ((P_PAave - PPUsyst_Nom1) / ...      % Nom pulm vasc resist
            CO_RHC) * (60/1000);                        %  with conv to mmHg*s/mL
    elseif (0.975 * P_PAsyst >= PPUsyst_Nom1)
        Rpul_Nom = ((P_PAsyst - PPUsyst_Nom1) / ...     % Nom pulm vasc resist
            CO_RHC) * (60/1000);                        %  with conv to mmHg*s/mL
    else
        Rpul_Nom = ((P_PAsyst - PPUsyst_Nom2) / ...     % Nom pulm vasc resist
            CO_RHC) * (60/1000);                        %  with conv to mmHg*s/mL
    end
    
    % Systemic arteries and systemic veins parameters
    Eessa_Nom = (P_AOsyst - P_AOdiast) / V_SAsyst;      % Nom systemic art elastance
    Eessv_Nom = PSVpp_Nom / V_SVsyst;                   % Nom systemic ven elstnce
    P_AOave = (P_AOsyst/3) + ((2*P_AOdiast)/3);
    if (P_RVdiast <= 0)                                 % Again if P_RVdiast is
        PSVdiast_Nom = 1.025 * PRVdiast_Nom;            %  negative in the data
        PSVsyst_Nom = PSVdiast_Nom + PSVpp_Nom;
    else
        PSVdiast_Nom = 1.025 * P_RVdiast;
        PSVsyst_Nom = PSVdiast_Nom + PSVpp_Nom;
    end
    
    Rsys_Nom = ((P_AOave - PSVsyst_Nom) / ...           % Nom systemic vasc resist
        CO_RHC) * (60/1000);                            %  with conv to mmHg*s/mL
    
    % Heart valve paramenters
    PPUdiast_Nom = P_PCWave - (PPUpp_Nom/3);            % Nom pulm venous dstlc P
    PLVdiast_Nom = 0.975 * PPUdiast_Nom;
    Rmt_Nom = ((PPUdiast_Nom - PLVdiast_Nom) / ...      % Nom mtrl vlv resist
        CO_RHC) * (60/1000);                            %  with conv to mmHg*s/mL
    Rav_Nom = ((PLVsyst_Nom - P_AOsyst) / ...           % Nom artc vlv resist
        CO_RHC) * (60/1000);                            %  with conv to mmHg*s/mL
    if (P_RVdiast <= 0)                                 % Again if P_RVdiast is
        Rtc_Nom = ((PSVdiast_Nom - PRVdiast_Nom) / ...  % Nom tricspd vlv resist
            CO_RHC) * (60/1000);                        %  with PRVdiast estimate
    else                                                %  and conv to mmHg*s/mL
        Rtc_Nom = ((PSVdiast_Nom - P_RVdiast) / ...     % Nom tricspd vlv resist
            CO_RHC) * (60/1000);                        %  using PRVdiast from data
    end                                                 %  and conv to mmHg*s/mL
    if (P_PAsyst >= P_RVsyst)                           % In some cases Ppasyst
        PPAsyst_Nom = 0.975 * P_RVsyst;                 %  is larger than Prvsyst
    else                                                %  which is not possible
        PPAsyst_Nom = P_PAsyst;                         %  so for calculation here
    end                                                 %  reset to 0.975 * Prvsyst
    Rpv_Nom = (((P_RVsyst - PPAsyst_Nom) / ...          % Nom Pulm vlv resist w/conv
        CO_RHC) * (60/1000)) / 5;                       %  and fact based on opt vals
    
    NomParam_Values = {Eeslvf_Nom P0lvf_Nom lambdalvf_Nom ...
        Eesrvf_Nom P0rvf_Nom lambdarvf_Nom Eespa_Nom Eespu_Nom ...
        Rpul_Nom Eessa_Nom Eessv_Nom Rsys_Nom Rmt_Nom Rav_Nom ...
        Rtc_Nom Rpv_Nom};
    NomParam_Fields = {'Eeslvf_Nom' 'P0lvf_Nom' 'lambdalvf_Nom' ...
        'Eesrvf_Nom' 'P0rvf_Nom' 'lambdarvf_Nom' 'Eespa_Nom' 'Eespu_Nom' ...
        'Rpul_Nom' 'Eessa_Nom' 'Eessv_Nom' 'Rsys_Nom' 'Rmt_Nom' 'Rav_Nom' ...
        'Rtc_Nom' 'Rpv_Nom'};
    NomParam_Struct = cell2struct(NomParam_Values,NomParam_Fields,2);
    
    
end

