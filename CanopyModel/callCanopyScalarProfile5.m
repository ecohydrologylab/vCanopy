%% Function to compute scalar profile within the canopy
function [Canopy] = callCanopyScalarProfile5(option,Weather, ...
    Canopy,Soil,Constants)

dz = diff(Canopy.z);
dz = [dz(1) dz];
dy = diff(Canopy.windProfile);
Km = Canopy.lm^2*abs(dy./dz);
Km = [Km(1) Km];
dKm_dz = diff(Km)./dz;
dKm_dz = [dKm_dz(1) dKm_dz];

totalPressure = (Weather.pressure + 0*Canopy.eaProfile); % [Pa]
molarDensity = @ (temperature,Pressure)  ((Pressure.*273.15)./(22.4*101.35*(temperature + 273.15))); % [mol m-3]
molarVolume = @ (temperature,Pressure)  ((22.4*101.35*(temperature + 273.15))./(Pressure.*273.15)); % [m3 mol-1]
latentHeatVapour = @(temperature) (2500 - 2.36*temperature)*18; % [J mol-1]
StrFactor = 1;

%% Stem-Sink
stemHeatCapacity = Canopy.drySM*(4.18E3/3 + 0.85/(1-0.85)*4.18E3);%/(Weather.deltaT*3600);

%% Initialize canopy wind profile [m s-1]
if strcmp(option,'temperature') % Temperature Profile
    cN = Weather.tAir; % [deg C]
    scalarProfileLag = Canopy.tAirProfileLag;
    if isnan(sum(scalarProfileLag(2:end-1))) 
        StrFactor = 0; scalarProfileLag(:) = 0; stemHeatCapacity(:) = 0;
    end
    B = Canopy.H(2:end)/Constants.Cp + ...
        stemHeatCapacity(2:end)/Constants.Cp.*scalarProfileLag(2:end-1)/(Weather.deltaT*3600); % Source/Sink strength [C s-1]
    F = Soil.H/Constants.Cp;
elseif strcmp(option,'ea') % Vapour Profile
    cN = Weather.ea/totalPressure(end); % [Pa] 
    scalarProfileLag = Canopy.airQProfileLag;
    stemHeatCapacity(:) = 0;
    if isnan(sum(scalarProfileLag(2:end-1)))
        StrFactor = 0; scalarProfileLag(:) = 0; stemHeatCapacity(:) = 0;
    end
    B = Canopy.LE(2:end)./latentHeatVapour(Canopy.tAirProfile(2:end-1));
    F = Soil.LE./latentHeatVapour(Canopy.tAirProfile(2));
elseif strcmp(option,'ca') % CO2 Profile
    cN = Weather.ca; % [ppm] or [umol mol-1]
    stemHeatCapacity(:) = 0; StrFactor = 0; scalarProfileLag(:) = 0;
    B = -Canopy.aNet(2:end); % [u moles mol-1 s-1] = [ppm s-1]
    F = Soil.CO2Flux;
end

scalarProfile = [NaN linspace(cN,cN,length(Canopy.z))];
B = [NaN B NaN];

%% Formulation of A Matrix and b vector of AX = B
aMatrix = zeros(length(Canopy.z));
bVector = zeros(length(Canopy.z),1);

%% Bottom Node: Flux boundary condition
aMatrix(1,1) = molarDensity(Canopy.tAirProfile(2),totalPressure(2))*...
    (StrFactor*dz(2)/(Weather.deltaT*3600) + Km(2)/dz(2)) + ...
    stemHeatCapacity(2)/(Weather.deltaT*3600*Constants.Cp);    
aMatrix(1,2) = -molarDensity(Canopy.tAirProfile(2),totalPressure(2))*Km(2)/dz(2);
bVector(1) = B(2) + (1 - dz(2)/Km(2)*dKm_dz(2))*F;


%% Interior Nodes
for zLoop=2:1:length(scalarProfile)-2
    aMatrix(zLoop,zLoop-1) = molarDensity(Canopy.tAirProfile(zLoop+1),totalPressure(zLoop+1))*...
        (dKm_dz(zLoop+1) - Km(zLoop+1)/dz(zLoop+1));
    aMatrix(zLoop,zLoop) = molarDensity(Canopy.tAirProfile(zLoop+1),totalPressure(zLoop+1))*...
        (- dKm_dz(zLoop+1) + 2*Km(zLoop+1)/dz(zLoop+1)) + ...
        stemHeatCapacity(zLoop+1)/(Weather.deltaT*3600*Constants.Cp);
    aMatrix(zLoop,zLoop+1) = -molarDensity(Canopy.tAirProfile(zLoop+1),totalPressure(zLoop+1))*Km(zLoop+1)/dz(zLoop+1);
    bVector(zLoop) = B(zLoop+1);
end
%% Top Node: Dirichlet boundary condition
aMatrix(end,end) = 1;
bVector(end) = cN;    

%% Computation
scalarProfile = [NaN (aMatrix\bVector)'];

%% Result
if strcmp(option,'temperature')
    Canopy.tAirProfile = scalarProfile; 
    Canopy.FHs = -molarDensity(Canopy.tAirProfile(end-1),totalPressure(end-1))*Km(end-1)*...
        Constants.Cp*(scalarProfile(end)-scalarProfile(end-1))/dz(end);   
elseif strcmp(option,'ea')
    Canopy.airQProfile = scalarProfile;
    Canopy.eaProfile = Canopy.airQProfile.*totalPressure;
    Canopy.FLE = -molarDensity(Canopy.tAirProfile(end-1),totalPressure(end-1))*Km(end-1)*...
         latentHeatVapour(Canopy.tAirProfile(end-1))*(scalarProfile(end)-scalarProfile(end-1))/dz(end);    
elseif strcmp(option,'ca')
    Canopy.caProfile = scalarProfile;
    Canopy.FCO2 = Km(end-1)*(scalarProfile(end)-scalarProfile(end-1))/dz(end)/...
        molarVolume(Canopy.tAirProfile(end-1),totalPressure(end-1));    
end

end
