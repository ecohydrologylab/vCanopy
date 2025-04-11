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

% Hourly canopy simulations simulations
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

% tLeaf = reshape([CanopyOut.tLeafProfile],[11,48]);
% tsunLeaf = reshape([CanopyOut.sunTleaf],[11,48]);
% tshadeLeaf = reshape([CanopyOut.shadeTleaf],[11,48]);

% figure(1)
% subplot(221)
% hold on
% surf(tLeaf)
% title("Tleaf")
% subplot(222)
% hold on
% surf(tsunLeaf)
% title("sun Tleaf")
% subplot(223)
% hold on
% surf(tshadeLeaf)
% title("shadeTleaf")
% subplot(224)
% hold on
% plot([CanopyOut.FHs])
% title("H")

% Storing daily canopy output
DailyOutput.CanopyOut = CanopyOut;
DailyOutput.logFile = logFile;
end