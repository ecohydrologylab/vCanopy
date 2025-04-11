function [Soil] = callSoilFlux(Constants,EnergyOptions,Canopy,Soil,Weather)

Constants.R0 = 1.4; % Soil respiration rate at 10^oC [mu mole m-2 s-1]
Constants.Q10 = 2; % Temperature sensitivity of soil respiration [-]

totalPressure = (Weather.pressure + Canopy.eaProfile);
callLv = @(temperature) (2500 - 2.36*temperature)*18;
dM = Constants.vonk^2/(log(Canopy.z(2)/Soil.roughnessHeight))^2; % Momentum transfer coefficient
dV = dM*0.622*Constants.airdensity/(totalPressure(2)*Constants.rhoW); % Vapor transfer coefficient
dH = dM*Constants.Cp*Constants.airdensity; % Heat transfer coefficient

callEs = @(temperature)...
    (0.611*exp(17.502*temperature/(240.97+temperature)))*1000; % Saturation Vapour Pressure [Pa]

% Soil net radiation flux
callNetRadiation = @(temperature) Soil.PARAbsorbed+Soil.NIRAbsorbed+ ...
    Soil.longAbsorbed- EnergyOptions.epsilonSoil*Constants.Boltzman*(temperature+273.15)^4;
Soil.netRadiation = callNetRadiation(Soil.temperature);

% Sensible heat flux
callSensibleHeat = @(temperature) dH*Canopy.windProfile(2)*...
    (temperature-Canopy.tAirProfile(2));
Soil.H = callSensibleHeat(Soil.temperature);

% Latent Heat Flux
RH_soil = Soil.RH/100; % Relative humidity at soil
callLatentHeat = @(temperature) dV*Canopy.windProfile(2)*...
    (callEs(temperature)*RH_soil - Canopy.eaProfile(2))* ...
    Constants.rhoW*callLv(temperature);
Soil.LE = callLatentHeat(Soil.temperature);

% Soil ground heat flux
callGFlux = @(temperature) (Constants.soilK *...
    (temperature-Canopy.tAirProfile(2))/(0.05));
Soil.G = callGFlux(Soil.temperature);

% Calculating soil temperature
fSoilEnergyFlux = @(temperature)(callNetRadiation(temperature) - ...
    (callSensibleHeat(temperature) + callLatentHeat(temperature) + ...
    callGFlux(temperature)))^2;
[temperature,fval,exitFalg] = fminbnd(fSoilEnergyFlux,4,50);

if exitFalg ~= 1
    disp("tSoil_junk") 
    temperature = Canopy.tAirProfile(end);
end
relax = 0.8;
Soil.temperature = temperature*relax + (1-relax)*Soil.temperature;

% Soil Sensible heat fluxe
Soil.H = callSensibleHeat(Soil.temperature);
% Soil Latent Heat Flux
Soil.LE = callLatentHeat(Soil.temperature);
% Soil ground heat flux
Soil.G = callGFlux(Soil.temperature);
% Soil CO2 Flux calculation
Soil.CO2Flux = Constants.R0* Constants.Q10^((Soil.temperature-10)/10);
% Soil net Radiation
Soil.netRadiation = callNetRadiation(Soil.temperature);
end