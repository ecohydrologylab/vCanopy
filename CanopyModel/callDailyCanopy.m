function [DailyOutput] = callDailyCanopy(Constants,EnergyOptions,LeafBoundaryLayer,...
    Photosynthesis,Stomata,Daily)

% Extracting the varaibles from Daily
Weather = Daily.Weather;
Canopy = Daily.Canopy;
Soil = Daily.Soil;
groupIndex = Daily.groupIndex;
logFile = string(Daily.logFile);

% Stored simualtions results in case of modified canopy simulations
callStoredSimulations

% Initializing canopy out
callInitializePreviousCanopy

% Applying season and vertical decay to photosynthesis
Photosynthesis = callModifyPhotosyntheticCapacity(Photosynthesis,...
    Weather(1),Canopy(1));

% Diurnal canopy simulations simulations
for hLoop = 1:1:length(Weather)
    
    % Compute canopy model
    [CanopyOut(hLoop),logFile(hLoop)] = callCanopy(...
        Constants,EnergyOptions,LeafBoundaryLayer,Photosynthesis,Stomata,...
        Weather(hLoop),Soil(hLoop),Canopy(hLoop),PreviousCanopy,savedtLeaf(hLoop,:),...
        savedcaProfile(hLoop,:),savedeaProfile(hLoop,:),...
        savedtAirProfile(hLoop,:),savedRHProfile(hLoop,:),...
        groupIndex(hLoop),savedlongAbsorbed(hLoop,:),savedlongEmitted(hLoop,:),...
        savedwindProfile(hLoop,:));
    
    % Updating the previous canopy
    callUpdatePreviousCanopy      
end

% Storing daily canopy output
DailyOutput.CanopyOut = CanopyOut;
DailyOutput.logFile = logFile;
end