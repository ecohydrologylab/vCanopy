%% Function to compute leaf temperature using energy balance approach
function [LeafState,LeafMassFlux,LeafEnergyFlux] = ...
    callEnergyBalance(Constants,CanopyLayer,EnergyOptions,Weather,LeafMassFlux,LeafState)

% Predefined functions
callLv = @(temperature) (2500 - 2.36*temperature)*18; % [J mol-1]
callEs = @(temperature) (0.611*exp(17.502*temperature/(240.97+temperature)))*1000; %% Clausius–Clapeyron relation Saturation Vapour Pressure [Pa]
callSensibleHeat = @(gbh,temperature) (EnergyOptions.Hfactor*Constants.Cp*gbh*(temperature-Weather.tAir)); % Sensible heat flux [W m-2]
callLatentHeat = @(gv,temperature) (EnergyOptions.LEfactor*callLv(temperature)* ...
    gv/(Weather.pressure)*(callEs(temperature)-Weather.ea)); % Latent heat flux [W m-2]
callEmission = @(temperature) (EnergyOptions.LWfactor*EnergyOptions.epsilonLeaf*Constants.Boltzman* ...
    (273.15+temperature)^4.0); % Long wave radiation flux emitted [W m-2]
callStorage = @(temperature) 0*(1.396E3 + 0.8/(1-0.8)*4.18E3)/(Weather.deltaT*3600)....
    *CanopyLayer.MPLA*(temperature - CanopyLayer.tLeafLag);
gbh = LeafState.gb*0.924;


if EnergyOptions.switch == 1 && EnergyOptions.longwaveSwitch == 0
    callEnergyBalanceResidual = @(temperature) (Weather.PAR + Weather.NIR + Weather.long...
        - 0.506*LeafMassFlux.aNet-callSensibleHeat(gbh,temperature)-callLatentHeat(LeafState.g,temperature)- ...
        callEmission(temperature) - callStorage(temperature))^2; % Compute energy balance residual [W m-2]
else 
    callEnergyBalanceResidual = @(temperature) (Weather.PAR + Weather.NIR + Weather.long...
        - 0.506*LeafMassFlux.aNet-callSensibleHeat(gbh,temperature)-callLatentHeat(LeafState.g,temperature)- ...
        (1-CanopyLayer.tauD)/CanopyLayer.deltaLAI*callEmission(temperature) - callStorage(temperature))^2; % Compute energy balance residual [W m-2]
end

if EnergyOptions.switch == 1
    [temperature,SSE,optimizationFlag] = fminbnd(callEnergyBalanceResidual,2,60); % Compute leaf temperature [Celsius]
else
    callEnergyBalanceResidual = @(temperature) NaN;
    temperature = CanopyLayer.savedtLeaf; % Turn off leaf energy balance
end

%% Update values
LeafState.tLeaf = temperature; % [Celsius]
LeafEnergyFlux.radiation = Weather.PAR + Weather.NIR + Weather.long; % Net radiation flux [W m-2]
LeafEnergyFlux.H = callSensibleHeat(gbh,temperature); % Sensible heat flux [W m-2]
LeafEnergyFlux.LE = callLatentHeat(LeafState.g,temperature); % Latent heat flux [W m-2]
LeafEnergyFlux.emission = callEmission(temperature); % Long wave radiation flux emitted [W m-2]
LeafEnergyFlux.storage = callStorage(temperature);
LeafEnergyFlux.residual = (callEnergyBalanceResidual(temperature))^0.5; % Net radiation flux [W m-2]
LeafMassFlux.transpiration = LeafEnergyFlux.LE/callLv(temperature)*1E6; % Leaf transpiration [u moles m-2 s-1]

end