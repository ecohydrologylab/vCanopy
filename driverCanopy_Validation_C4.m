clc
clear variables
% close all
addpath('./CanopyModel')
addpath('./LeafModel')

%% Input
GrowthData = readtable("./Input/driverDrewryInput_C4_New.xlsx");
hour = unique(GrowthData.hour);     
%% SIP reduction
tic
Reduction = 0;

%% Computation
for index = 1:length(Reduction)
    
    %% Getting the Input
    canopyLayers = 20;
    Options.canopyPoints = canopyLayers+1; % Including the soil layer
    Options.plant = 4; % For C4 crops
    
    %% Light attenuation controls
    Options.isotropy = "UOC"; % Uniform overcast sky
    Options.diffuseFactor = 0; % No external influence on the diffusive fraction
    
    %% Switches
    Options.uniformLAD = 0; % If the leaf area index within the canopy
    Options.LEBSwitch = 1; % Energy balance switch in each layer
    Options.mixingSwitch = 1; % MicroEnvironment On: 1; Off: 0
    Options.nitrogenDecay = 1; % Nitrogen decay switch
    Options.longwaveSwitch = 1; % Whether longwave radiation is modeled
    Options.fixedMicro = 0; % This will fix the microenvironment/longwave but allow the leaf energu balance to work
    
    %% Changes based on crop species
    Options.GDDpre = 660; % Heat requirements during emergence to flowering
    Options.GDDpost = 850; % Heat requirements during flowering to maturity (Absolute = 1510)
    Options.X = 1.37; % Leaf angle distribution
    Options.Jmax_Vcmax = 5; % Jmax/Vcmax 
    Options.kn = 0.5; % Exponential nitrogen decay
    Options.PcSeasonalVariation = 1; % Decay in photosynthetic parameters
    Options.linearDTD = 0.471; % daily decline in Vcmax after 600 GDD
    Options.omega = 1; % Clumping factor for maize
    Options.aPAR = 0.80; % absorbtivity of PAR
    Options.aNIR = 0.23; % absorbtivity of NIR
    Options.epsilonLeaf = 0.94; % emissivity of leaf for LW
    
    %% Global Change simulations
    Options.LAIRatio = 1; % Fraction to model LAI increase under elevated [CO2]
    Options.vcmaxRatio = 1; % Factor to model photosynthesis capacity depreciation under elevated [CO2]
    Options.deltaTair = 0; % 2.7; % Increase in Tair for future climate
    Options.deltaRH = 0; % -3.5; % Percentage decrease in RH for future climate
    Options.deltaCO2 = 0; % 130; % Increase in [CO2] for future climate
    
    %% Input and initialization of canopy variables
    [Constants,EnergyOptions,LeafBoundaryLayer,Weather,Soil,Canopy,...
        Photosynthesis,Stomata] = callInputCanopy(Options,GrowthData);
    
    %% Sorting the data into days
    for dayLoop = 1:1:length(Weather)/length(hour)
        groupIndex = (dayLoop-1)*length(hour)+1:dayLoop*length(hour);
        Daily(dayLoop).Weather = Weather(groupIndex);
        Daily(dayLoop).Canopy = Canopy(groupIndex);
        Daily(dayLoop).Soil = Soil(groupIndex);
        Daily(dayLoop).groupIndex = groupIndex';
        Daily(dayLoop).logFile = groupIndex*NaN;
    end
    
    %% Applying stomatal reduction
    Stomata.slope = Stomata.slope*(1-Reduction(index)/100);
    Stomata.intercept = Stomata.intercept*(1-Reduction(index)/100);
    
    %% Setup output and log file name
    fileLocation = './Output/';
    fileName = strcat("val_vCanopy1D_","C",...
        num2str(Photosynthesis.plant(1)),"_soilRH_",num2str(Soil(1).RH),...
        "%_Ca_",num2str(Weather(1).ca),"_Vc_",num2str(Photosynthesis.vcmax25(end))...
        ,"_JxE1_",num2str(Photosynthesis.jmax25(end)*10));
    file = fopen(strcat(fileLocation,"LogFile_",fileName,'.txt'),'wt');
    
    %% Model Calculations
    parfor dayLoop = 1:height(GrowthData)/length(hour)
        DailyOutput(dayLoop) = callDailyCanopy(Constants,EnergyOptions,LeafBoundaryLayer,...
            Photosynthesis,Stomata,Daily(dayLoop));
    end
    
    %% Sorting back the daily output into table
    for dayLoop = 1:1:length(Weather)/length(hour)
        groupIndex = (dayLoop-1)*length(hour)+1:dayLoop*length(hour);
        CanopyOut(groupIndex,:) = struct2table(DailyOutput(dayLoop).CanopyOut);
        logFile(groupIndex) = DailyOutput(dayLoop).logFile;
    end
    
    clear Daily DailyOutput
    
    % Write and Save the log file
    for fLoop = 1:height(GrowthData)
        fprintf(file,logFile(fLoop));
    end
    fclose(file);
    
    %% Save the results
    save(strcat(fileLocation,fileName,".mat"));
end
toc
