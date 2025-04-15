function [Canopy] = callCanopyMicroenvironmentConstant(Constants,EnergyOptions,Weather,Canopy,Error,savedcaProfile,savedeaProfile,...
    savedtAirProfile,savedRHProfile,CanopyIteration)

if EnergyOptions.fixedMicro == 1
    Canopy.caProfile = savedcaProfile;
    Canopy.eaProfile = savedeaProfile;
    Canopy.tAirProfile = savedtAirProfile;
    Canopy.RHProfile = savedRHProfile;
else
    Canopy.tAirProfile = Weather.tAir*ones(1,length(Canopy.z)+1);
    Canopy.eaProfile = Weather.ea*ones(1,length(Canopy.z)+1);
    Canopy.caProfile = Weather.ca*ones(1,length(Canopy.z)+1);
    Canopy = callCanopyRHProfile(Canopy); % Compute RH profile
end

dz = diff(Canopy.z);
dz = [dz(1) dz];
dy = diff(Canopy.windProfile);
Km = Canopy.lm^2*abs(dy./dz);
Km = [Km(1) Km];

totalPressure = (Weather.pressure + Canopy.eaProfile);
molarVolume = @ (temperature,Pressure)  ((22.4*101.35*(temperature + 273.15))./(Pressure.*273.15));
latentHeatVapour = @(temperature) (2500 - 2.36*temperature)*18;

Canopy.FHs = -Km(end-1)*(Canopy.tAirProfile(end)-Canopy.tAirProfile(end-1))/dz(end)/...
    molarVolume(Canopy.tAirProfile(end-1),totalPressure(end-1))*Constants.Cp;
Canopy.FLE = -Km(end-1)*(Canopy.eaProfile(end)-Canopy.eaProfile(end-1))/dz(end)/...
    molarVolume(Canopy.tAirProfile(end-1),totalPressure(end-1))...
    *latentHeatVapour(Canopy.tAirProfile(end-1))/totalPressure(end);
Canopy.FCO2 = Km(end-1)*(Canopy.caProfile(end)-Canopy.caProfile(end-1))/dz(end)/...
    molarVolume(Canopy.tAirProfile(end-1),totalPressure(end-1));

if EnergyOptions.switch == 1
    Canopy.sunTleaf(1:length(Canopy.z)) =  Error.relax*...
        Canopy.sunTleaf+ (1-Error.relax)*CanopyIteration.sunTleaf;
    Canopy.shadeTleaf(1:length(Canopy.z)) = (Error.relax*...
        Canopy.shadeTleaf + (1-Error.relax)*CanopyIteration.shadeTleaf);
end

end