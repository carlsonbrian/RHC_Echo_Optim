% ***********************************************************************************
%                S I X   P A N E L   F I G U R E   F U N C T I O N   
% ***********************************************************************************
%
%   This function plots out a six panel figure showing 
%
%   Model originally created on     11 December 2018
%   Model last modfied on           25     July 2019
%
%   Reproduced by       Brian Carlson
%                       Physiological Systems Dynamics Laboratory
%                       Department of Molecular and Integrative Physiology
%                       University of Michigan
%
% ***********************************************************************************
%  START OF         S I X   P A N E L   F I G U R E   F U N C T I O N   
% ***********************************************************************************

    function SixPanel_Figure(AllStruct_Struct,SimPlot_Struct)
    
    %% UNPACKING THE STRUCTURES
    FlagData_Struct = AllStruct_Struct.FlagData_Struct;
    % Unpacking the option flags and patient data if present
    SmithParam_Flag = FlagData_Struct.SmithParam_Flag;
    Comp2Smith_Flag = FlagData_Struct.Comp2Smith_Flag;
    RHCEcho_Flag = FlagData_Struct.RHCEcho_Flag;
    PatData_Struct = AllStruct_Struct.PatData_Struct;
    % If we have RHC or RHC and Echo data extract those structures and values
    if (SmithParam_Flag == 0)
        RHCData_Struct = ...
            AllStruct_Struct.RHCData_Struct;
        DID_Num = PatData_Struct.DID_Num;
        P_RVsyst = RHCData_Struct.P_RVsyst;             % Systol RV pressure (mmHg) 
        P_RVdiast = RHCData_Struct.P_RVdiast;           % Diast RV pressure (mmHg)
        P_PAsyst = RHCData_Struct.P_PAsyst;             % Syst pulm art press (mmHg)
        P_PAdiast = RHCData_Struct.P_PAdiast;           % Diast pulm art press (mmHg)
        P_PCWave = RHCData_Struct.P_PCWave;             % Ave pulm wedge press (mmHg)
        P_AOsyst = RHCData_Struct.P_AOsyst;             % Systol aortic press (mmHg)
        P_AOdiast = RHCData_Struct.P_AOdiast;           % Diast aortic press (mmHg)
        CO_Fick = RHCData_Struct.CO_Fick;               % Crdc outpt Fick (L/min)
        CO_Thermo = RHCData_Struct.CO_Thermo;           % Crdc outpt thermo (L/min)
        if (RHCEcho_Flag == 1)
            EchoData_Struct = ...
                AllStruct_Struct.EchoData_Struct;
            ID_LVsyst = EchoData_Struct.ID_LVsyst;      % Systol LV inner diam (mm)
            ID_LVdiast = EchoData_Struct.ID_LVdiast;    % Diast LV inner diam (mm)
            HR_Echo = EchoData_Struct.HR_Echo;
            CO_EchoD = EchoData_Struct.CO_EchoD;        % Cardiac output Echo-Dop (L/min)
            V_LVsyst = EchoData_Struct.V_LVsyst;        % Systolic LV volume (mL)
            V_LVdiast = EchoData_Struct.V_LVdiast;      % Diastolic LV volume (mL)
            if (V_LVsyst == -1)
                V_LVsyst_Teich = ...                    % Calculate the 
                    ((ID_LVsyst/10)^3) / ...            %  LV systolic volume
                    ((6/pi) * ((0.075 * ...             %  using the Teichholz
                    (ID_LVsyst/10)) + 0.18));           %  formula
                V_LVdiast_Teich = ...                   % Calculate the 
                    ((ID_LVdiast/10)^3) / ...           %  LV diastolic volume     
                    ((6/pi) * ((0.075 * ...             %  using the Teichholz    
                    (ID_LVdiast/10)) + 0.18));          %  formula
                CO_EchoTeich = ((V_LVdiast_Teich - ...  % Calculate SV * HR
                    V_LVsyst_Teich) * ...               %  estimate of cardiac
                    HR_Echo) / 1000;                    %  output for Teich vols
            else
                CO_EchoSVHR = ((V_LVdiast - ...
                    V_LVsyst) * HR_Echo) / 1000;        % Crdc output SV * HR (L/min)
            end
        end
    end
    % Load the previously saved Smith data simulation
    %  when we want to compare patient data/sim to Smith
    if (Comp2Smith_Flag == 1)
        load('SmithVars.mat','T_SimSmith','P_RVSimSmith', ...
            'P_AOSimSmith','P_PASimSmith','P_PUSimSmith', ...
            'CO_SimSmith','P_LVSimSmith','P_VCSimSmith', ...
            'V_LVSimSmith','V_RVSimSmith');
    end
    % Unpacking the simulation results
    if (SmithParam_Flag == 1)
 
        T_Out = SimPlot_Struct.T_Out;
        P_RVSim = SimPlot_Struct.P_RVSim;
        P_AOSim = SimPlot_Struct.P_AOSim;
        P_PASim = SimPlot_Struct.P_PASim;
        P_PUSim = SimPlot_Struct.P_PUSim;
        CO_Sim = SimPlot_Struct.CO_Sim;
        V_LVSim = SimPlot_Struct.V_LVEchoSim;
        V_RVSim = SimPlot_Struct.V_RVEchoSim;
        P_LVSim = SimPlot_Struct.P_LVSim;
        P_VCSim = SimPlot_Struct.P_VCSim;
        
    else
        
        if (RHCEcho_Flag == 1)
            
            T_Out = SimPlot_Struct.T_OutRHC;
            T_OutEcho = SimPlot_Struct.T_OutEcho;
            P_RVSim = SimPlot_Struct.P_RVRHCSim;
            P_RVEchoSim = SimPlot_Struct.P_RVEchoSim;
            P_AOSim = SimPlot_Struct.P_AOSim;
            P_PASim = SimPlot_Struct.P_PASim;
            P_PUSim = SimPlot_Struct.P_PUSim;
            CO_RHCSim = SimPlot_Struct.CO_RHCSim;
            CO_EchoSim = SimPlot_Struct.CO_EchoSim;
            V_LVSim = SimPlot_Struct.V_LVEchoSim;
            V_RVSim = SimPlot_Struct.V_RVEchoSim;
            P_LVSim = SimPlot_Struct.P_LVRHCSim;
            P_LVEchoSim = SimPlot_Struct.P_LVEchoSim;
            P_VCSim = SimPlot_Struct.P_VCSim;
            
        else
            
            T_Out = SimPlot_Struct.T_Out;
            P_RVSim = SimPlot_Struct.P_RVSim;
            P_AOSim = SimPlot_Struct.P_AOSim;
            P_PASim = SimPlot_Struct.P_PASim;
            P_PUSim = SimPlot_Struct.P_PUSim;
            CO_RHCSim = SimPlot_Struct.CO_RHCSim;
            V_LVSim = SimPlot_Struct.V_LVSim;
            V_RVSim = SimPlot_Struct.V_RVSim;
            P_LVSim = SimPlot_Struct.P_LVSim;
            P_VCSim = SimPlot_Struct.P_VCSim;
            
        end
        
    end
    
    %% BUILDING THE FIGURE
    ScrSize = get(0,'ScreenSize');                      % Getting screen size
    
    SixP_Fig = figure('Position', ...                   % Positioning the figure
        [ScrSize(3)/50 ScrSize(4)/50 ...                %  on the screen
        ScrSize(3)/1.05 ScrSize(4)/1.05]); 
    % Figuring out which plot title to display
    if (SmithParam_Flag == 1)
        SupT_Hndl = suptitle({'Smith Params Normal'});
    elseif (Comp2Smith_Flag == 1)
        SupT_Hndl = suptitle({'Comparison to Smith Params'});
    else
        SupT_Hndl = suptitle({['DID Number ' DID_Num]});
    end
    set(SupT_Hndl,'FontSize',24,'FontWeight','bold')
    

    %% SUBPLOT OF P_RV, P_AO, P_PA AND P_PU AND COMP TO DATA IF PRESENT
    subplot(2,3,1)
    
    plot(T_Out,P_RVSim,'-g', ...                        % P_RV simulation
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
        plot([T_Out(1),T_Out(end)], ...                 % P_RV systole data
            [P_RVsyst,P_RVsyst],'-.g', ...
            'LineWidth',1.5, ...
            'HandleVisibility','off')
        plot([T_Out(1),T_Out(end)], ...                 % P_RV diastole data
            [P_RVdiast,P_RVdiast],':g', ...
            'LineWidth',1.5, ...
            'HandleVisibility','off')
    end
    plot(T_Out,P_AOSim,'-r', ...                        % P_AO simulation
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
        plot([T_Out(1),T_Out(end)], ...                 % P_AO systole data
            [P_AOsyst,P_AOsyst],'-.r', ...
            'LineWidth',1.5, ...
            'HandleVisibility','off')
        plot([T_Out(1),T_Out(end)], ...                 % P_AO diastole data
            [P_AOdiast,P_AOdiast],':r', ...
            'LineWidth',1.5, ...
            'HandleVisibility','off')
    end
    plot(T_Out,P_PASim,'-b', ...                        % P_PA simulation
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
        plot([T_Out(1),T_Out(end)], ...                 % P_PA systole data
            [P_PAsyst,P_PAsyst],'-.b', ...
            'LineWidth',1.5, ...
            'HandleVisibility','off')
        plot([T_Out(1),T_Out(end)], ...                 % P_PA diastole data
            [P_PAdiast,P_PAdiast],':b', ...
            'LineWidth',1.5, ...
            'HandleVisibility','off')
    end
    plot(T_Out,P_PUSim,'-c', ...                        % P_PU simulation
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
        plot([T_Out(1),T_Out(end)], ...                 % P_PCW average data
            [P_PCWave,P_PCWave],'-.c', ...
            'LineWidth',1.5, ...
            'HandleVisibility','off')
    end
    if (SmithParam_Flag == 0)
        PMax_SP1 = 1.45 * ...
            max(P_AOsyst,max(P_AOSim));
        if (CO_Thermo ~= -1)
            text(0.5,0.95*PMax_SP1, ...                 % CO RHC Thermo data
                ['CO RHC Data = ' num2str(CO_Thermo)])
        else
            text(0.5,0.95*PMax_SP1, ...                 % CO RHC Fick data
                ['CO RHC Data = ' num2str(CO_Fick)])
        end
        text(0.5,0.90*PMax_SP1, ...                     % CO RHC simulation
            ['CO RHC Sim  = ' num2str(CO_RHCSim)]) 
        if (RHCEcho_Flag == 1)            
            if (V_LVsyst == -1)
                text(0.5,0.85*PMax_SP1, ...             % CO Echo Teich data
                    ['CO Teich Data = ' ...
                    num2str(CO_EchoTeich)])
            else
                text(0.5,0.85*PMax_SP1, ...             % CO Echo SVHR data
                    ['CO SVHR Data = ' ...
                    num2str(CO_EchoSVHR)])
            end
            if (CO_EchoD ~= -1)
                text(0.5,0.80*PMax_SP1, ...             % CO EchoD data
                    ['CO EchoD Data = ' ...
                    num2str(CO_EchoD)])
                text(0.5,0.75*PMax_SP1, ...             % CO Echo simulation
                    ['CO Echo Sim  = ' ...
                    num2str(CO_EchoSim)])
            else
                text(0.5,0.80*PMax_SP1, ...             % CO Echo simulation
                    ['CO Echo Sim  = ' ...
                    num2str(CO_EchoSim)])
            end
        end
        
    elseif (Comp2Smith_Flag == 1)
        PMax_SP1 = 1.45 * ...
            max([max(P_AOSim) max(P_AOSimSmith)]);
        text(0.5,0.95*PMax_SP1, ...                     % CO Smith params
            ['CO Smith = ' num2str(CO_SimSmith)])
        text(0.5,0.90*PMax_SP1, ...
            ['CO Sim  = ' num2str(CO_Sim)])             % CO simulation
    else
        PMax_SP1 = 1.45 * max(P_AOSim);
        text(0.5,0.90*PMax_SP1, ...
            ['CO Sim  = ' num2str(CO_Sim)])             % CO simulation
    end
    xlim([0 5])                                         % Formatting subplot
    ylim([-20 PMax_SP1])
    LegHndl1 = legend('show');
    set(LegHndl1,'Box','off','FontSize',8)
    set(gca,'FontSize',14,'FontWeight', ...
        'bold','Box','off')
    xlabel('Time (sec)','FontSize',20, ...
        'FontWeight','bold')
    ylabel('Pressures (mmHg)', ...
        'FontSize',20,'FontWeight','bold')
    

    %% SUBPLOT OF V_LV and V_RV WITH V_LV DATA SHOWN
    subplot(2,3,2)
        if (RHCEcho_Flag == 1 && SmithParam_Flag == 0)
            T_OutVols = T_OutEcho;
        else
            T_OutVols = T_Out;
        end
        plot(T_OutVols,V_LVSim,'-k', ...                % V_LV simulation    
            'LineWidth',3)
        hold on
        if (Comp2Smith_Flag == 1)
            SmithP5 = plot(T_SimSmith, ...
                V_LVSimSmith,'--k', ...
                'LineWidth',3, ...
                'HandleVisibility','off');
            SmithP5.Color(4) = 0.25;
        end
        plot(T_OutVols,V_RVSim,'-g', ...                % V_RV simulation 
            'LineWidth',3)
        if (Comp2Smith_Flag == 1)
            SmithP6 = plot(T_SimSmith, ...
                V_RVSimSmith,'--g', ...
                'LineWidth',3, ...
                'HandleVisibility','off');
            SmithP6.Color(4) = 0.25;
        end
        if (Comp2Smith_Flag == 1)
            VMax_SP2 = 1.20*max([max(V_LVSim) ...
                max(V_RVSim) ...
                max(V_LVSimSmith)]);
            VMin_SP2 = 0.85*min([min(V_LVSim) ...
                min(V_RVSim) ...
                min(V_LVSimSmith)]);
        else
            VMax_SP2 = 1.20*max(max(V_LVSim), ...
                max(V_RVSim));
            VMin_SP2 = 0.85*min(min(V_LVSim), ...
                min(V_RVSim));
        end
        if ((SmithParam_Flag == 0) && (RHCEcho_Flag == 1))
            if (V_LVsyst ~= -1)
                T_OutVols(1)
                T_OutVols(end)
                plot([T_OutVols(1) T_OutVols(end)], ...
                    [V_LVsyst V_LVsyst], ...
                    '-.k','LineWidth',1.5)
                plot([T_OutVols(1) T_OutVols(end)], ...
                    [V_LVdiast V_LVdiast], ...
                    ':k','LineWidth',1.5)
            else
                plot([T_OutVols(1) T_OutVols(end)], ...
                    [V_LVsyst_Teich V_LVsyst_Teich], ...
                    '-.k','LineWidth',1.5)
                plot([T_OutVols(1) T_OutVols(end)], ...
                    [V_LVdiast_Teich V_LVdiast_Teich], ...
                    ':k','LineWidth',1.5)
            end
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
        

    %% SUBPLOT OF P_LV, P_AO, P_PU and P_VC WITHOUT DATA
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
        

    %% SUBPLOT OF P_RV, P_PA and P_VC WITH NO DATA
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
        

    %% SUBPLOT OF LV AND RV PV LOOPS WITH DATA IF ECHO DATA AVAILABLE  
    subplot(2,3,[3 6]) 
        LVVol_Max = max(V_LVSim);
        RVVol_Max = max(V_RVSim);
        if (SmithParam_Flag == 0)
            if (Comp2Smith_Flag == 1)
                LVVol_SmithMax = max(V_LVSimSmith);
                if (RHCEcho_Flag == 1)
                    VMax_SP36 = max([LVVol_Max RVVol_Max ...
                        V_LVdiast LVVol_SmithMax]);
                else
                    VMax_SP36 = max([LVVol_Max RVVol_Max ...
                        LVVol_SmithMax]);
                end
            else
                if (RHCEcho_Flag == 1)
                    VMax_SP36 = max([LVVol_Max RVVol_Max V_LVdiast]);
                else
                    VMax_SP36 = max([LVVol_Max RVVol_Max]);
                end
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
        if (RHCEcho_Flag == 1 && SmithParam_Flag == 0)
            P_LVSim = P_LVEchoSim;
        end
        plot(V_LVSim,P_LVSim,'-k', ...              % V_LV vs P_LV simulation
            'LineWidth',3)
        hold on
        if (Comp2Smith_Flag == 1)
            SmithP13 = plot(V_LVSimSmith, ...
                P_LVSimSmith,'--k', ...
                'LineWidth',3, ...
                'HandleVisibility','off');
            SmithP13.Color(4) = 0.15;
        end
        if (RHCEcho_Flag == 1 && SmithParam_Flag == 0)
            P_RVSim = P_RVEchoSim;
        end
        plot(V_RVSim,P_RVSim,'-g', ...              % V_RV vs P_RV simulation
            'LineWidth',3)
        if (Comp2Smith_Flag == 1)
            SmithP14 = plot(V_RVSimSmith, ...
                P_RVSimSmith,'--g', ...
                'LineWidth',3, ...
                'HandleVisibility','off');
            SmithP14.Color(4) = 0.15;
        end
        if ((SmithParam_Flag == 0) && (RHCEcho_Flag == 1))
            if (V_LVsyst ~= -1)
                plot([V_LVsyst V_LVsyst], ...
                    [0.5*PMin_SP36 1.05*PMax_SP36], ...
                    '-.k','LineWidth',1.5)
                plot([V_LVdiast V_LVdiast], ...
                    [0.5*PMin_SP36 1.05*PMax_SP36], ...
                    ':k','LineWidth',1.5)
            else
                plot([V_LVsyst_Teich V_LVsyst_Teich], ...
                    [0.5*PMin_SP36 1.05*PMax_SP36], ...
                    '-.k','LineWidth',1.5)
                plot([V_LVdiast_Teich V_LVdiast_Teich], ...
                    [0.5*PMin_SP36 1.05*PMax_SP36], ...
                    ':k','LineWidth',1.5)
            end
        end
        xlim([0 VMax_SP36*1.20])                % Formatting subplot
        ylim([PMin_SP36 1.15*PMax_SP36])
        LegHndl36 = legend('LV','RV');
        set(LegHndl36,'Box','off','FontSize',8)
        set(gca,'FontSize',14, ...
            'FontWeight','bold','Box','off')
        xlabel('LV/RV Volume (mL)', ...
            'FontSize',20,'FontWeight','bold')
        ylabel('LV/RV Pressure (mmHg)', ...
            'FontSize',20,'FontWeight','bold')
       
        % Saving figure as a png 
        print(SixP_Fig,'-dpng','SixPanel_Fig.png')
%         hgsave(SixP_Fig,'SixPanel_Fig.png')
        
end

