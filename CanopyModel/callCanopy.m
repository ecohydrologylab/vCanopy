function [CanopyOut,logFile] = ...
    callCanopy(Constants,EnergyOptions,LeafBoundaryLayer,Photosynthesis,...
    Stomata,Weather,Soil,Canopy,PreviousCanopy,savedtLeaf,savedcaProfile,savedeaProfile,...
    savedtAirProfile,savedRHProfile,loop,savedlongAbsorbed,...
    savedlongEmitted,savedwindProfile)

addpath('./CanopyModel')
addpath('./LeafModel')
nanFlag = 1;

% Compute sun and shade fractions
Canopy = callSunShadeFraction(Weather,Canopy);

% Compute short radiation attenuation
[Soil,Canopy] = callShortRadiationAttenuation(EnergyOptions,Weather,Soil,Canopy);

% Compute wind profile
nLayer = length(Canopy.z);
if EnergyOptions.mixingSwitch == 1 && nLayer > 2
    Canopy = callCanopyWindProfile(Canopy,Weather,Soil); % Compute wind profile
elseif EnergyOptions.mixingSwitch == 0 && EnergyOptions.fixedMicro == 1
    Canopy.windProfile = savedwindProfile;
else
    Canopy.windProfile = zeros(1,size(Canopy.z,2)+1)+Weather.wind;
end

% Assign and initialize canopy and leaf fluxes and states
callBetterInitializeCanopy
[Photosynthesis,LeafState,LeafMassFlux,LeafEnergyFlux] = ...
    callOneTimeLeafInitialization(Photosynthesis,Weather,Initialize);

% Convergence loop for longwave and leaf temperature
Error.maxLoop = 50;
Error.tolerance = 0.1;
Error.sunTleaf = 10;
Error.shadeTleaf = 10;
Error.caProfile = 10;
Error.eaProfile = 10;
Error.tAirProfile = 10;
Error.soilTemperature = 10;
Error.relax = 0.5;
Error.loop = 1;

while (Error.loop < Error.maxLoop) && ...
        ((Error.sunTleaf > Error.tolerance) || ...
        (Error.shadeTleaf > Error.tolerance) || ...        
        (Error.caProfile > Error.tolerance) || ...
        (Error.eaProfile > Error.tolerance) || ...
        (Error.tAirProfile > Error.tolerance) || ...
        (Error.soilTemperature > Error.tolerance))
    
    % Compute long radiation attenuation
    if EnergyOptions.longwaveSwitch == 1
        [Soil,Canopy] = callLongRadiationAttenuation(Constants,EnergyOptions,Weather,Soil,Canopy);
    elseif (EnergyOptions.switch == 0 || EnergyOptions.fixedMicro == 1) && EnergyOptions.longwaveSwitch == 0
        Canopy.tauD = exp(-Canopy.kD.*Canopy.omega.*Canopy.deltaLAI);
        Canopy.longEmitted = savedlongEmitted;
        Canopy.longAbsorbed = savedlongAbsorbed;
        Canopy.longUp = savedlongAbsorbed*NaN;
        Canopy.sunLongAbsorbed = Canopy.sunFraction.*Canopy.longAbsorbed;
        Canopy.shadeLongAbsorbed = Canopy.shadeFraction.*Canopy.longAbsorbed;
        Soil.longAbsorbed = Canopy.longAbsorbed(1);
        Soil.longEmitted = Canopy.longEmitted(1);
    end
    % In case modified longwave calculations are required
    if EnergyOptions.longwaveSwitch == 0 && EnergyOptions.switch == 1 && EnergyOptions.fixedMicro == 0
        [Soil,Canopy] = callLongRadiationConstant(Constants,EnergyOptions,Weather,Soil,Canopy);
    end    

    % Assemble input for leaf model
    [LeafWeather,CanopyLayer] = callLeafInput(Weather,Canopy,savedtLeaf);
    
    % Storing old canopy/soil for canopy iteration
    CanopyIteration = Canopy;
    SoilTemp = Soil.temperature;
    
    % Simulating layers within canopy
    for leafLoop = 1:2*length(Canopy.z)
        if (LeafWeather(leafLoop).fraction ~= 0)
            [LeafMassFlux(leafLoop),LeafEnergyFlux(leafLoop),...
                LeafState(leafLoop)] = callNewtonLeaf4(Constants,...
                LeafBoundaryLayer,EnergyOptions,Photosynthesis(leafLoop),Stomata,...
                LeafWeather(leafLoop),CanopyLayer(leafLoop),LeafState(leafLoop),...
                LeafMassFlux(leafLoop),LeafEnergyFlux(leafLoop),loop);
        end
    end    
    
    Canopy.sunTleaf = [LeafState(1:nLayer).tLeaf];
    Canopy.shadeTleaf = [LeafState(nLayer+1:end).tLeaf];
    Canopy.tLeafProfile = Canopy.sunTleaf.*Canopy.sunFraction + Canopy.shadeTleaf.*Canopy.shadeFraction;
    
    % Compute canopy source/sink fluxes from leaf mass and energy fluxes
    [Canopy] = callCanopySourceSink(LeafMassFlux,LeafEnergyFlux,Canopy,Soil);    
    
    % Soil flux calculations
    [Soil] = callSoilFlux(Constants,EnergyOptions,Canopy,Soil,Weather);

    % Canopy microenvironment simulations (Turbulent mixing)
    if EnergyOptions.mixingSwitch == 1 && nLayer > 2
        [Canopy] = callCanopyMicroenvironment(EnergyOptions,Constants,Weather,Soil...
            ,Canopy,Error,CanopyIteration);
    else
        [Canopy] = callCanopyMicroenvironmentConstant(Constants,EnergyOptions,Weather,Canopy,Error,...
            savedcaProfile,savedeaProfile,savedtAirProfile,...
            savedRHProfile,CanopyIteration);
    end
        
    % Convergence Check
    [Canopy,Error] = callComputeError(Error,Canopy,CanopyIteration);
    Error.soilTemperature = abs(Soil.temperature-SoilTemp);
       
    % This wil check if its worth doing further canopy iterations
    if ~isreal(sum(Canopy.eaProfile(2:end))) || ~isreal(sum(Canopy.tAirProfile(2:end))) ...
        || ~isreal(sum(Canopy.caProfile(2:end))) || (Error.loop > 5 && sum(Canopy.RHProfile > 110) > 5) ...
        || sum(isnan(Canopy.caProfile(2:end))) || sum(isnan(Canopy.eaProfile(2:end))) ...
        || sum(isnan(Canopy.tAirProfile(2:end)))
        nanFlag = NaN;
        Error.loop = 1000;
    end
   
end

% Assembling final canopy output
callFinalOutput

end