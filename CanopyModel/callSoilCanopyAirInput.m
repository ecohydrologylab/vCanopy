%% Function to create canopy input
function [Soil,Canopy] = callSoilCanopyAirInput(Weather,Options)

%% Soil properties
Soil = table;
Soil.temperature = Weather.tAir+1; % Soil temperature [Celsius]
Soil.wind = zeros(height(Weather),1); % Wind speed at soil [m s-1]
Soil.RH = 85+zeros(height(Weather),1); % Relative humidity [%]
Soil.roughnessHeight = 0.005*ones(height(Weather),1); % Soil sroughness height [m]
Soil.CO2Flux = 1*ones(height(Weather),1); % Soil respiration [u mol m-2 s-1]
Soil.LE = 1*ones(height(Weather),1); % Soil latent heat flux [W m-2]
Soil.H = 1*ones(height(Weather),1); % Soil sensible heat flux [W m-2]
Soil.G = 1*ones(height(Weather),1); % Ground heat flux [W m-2]
Soil.netRadiation = zeros(height(Weather),1); % Net radiation absorbed by soil [W m-2]

%% Canopy properties
Canopy = table;
try
    Canopy.hN = Weather.Height; % Canopy height [m]
catch
    Canopy.hN = 2.5*ones(height(Weather),1); % Canopy height [m]
end

Canopy.ID = [1:1:height(Weather)]';
deltaZ = Canopy.hN/(Options.canopyPoints-1);
Canopy.deltaZ = [-1*deltaZ/2,repmat(deltaZ,[1,Options.canopyPoints-1])];
Canopy.z = cumsum(Canopy.deltaZ,2);
Canopy.deltaZ(:,1) = NaN; % Canopy vertical descritization unit [m]

Canopy.lm = 0.1367 * Canopy.hN; % Prandtl mixing-length [m]
Canopy.cd = 0.2*ones(height(Weather),1); % Canopy drag coefficient [-]
Canopy.x = Options.X*ones(height(Weather),1); % Leaf angle distribution factor [-]
Canopy.omega = Options.omega*ones(height(Weather),1); % Clumping factor [-]

if Options.canopyPoints == 2 % For big leaf model
    Canopy.deltaLAI = [Weather.LAI*0 Options.LAIRatio*Weather.LAI]; % Add soil layer
    Canopy.LAI = Canopy.deltaLAI;
else
    if Options.uniformLAD == 1
        LADnorm = ones(size(Canopy.z,1),size(Canopy.z,2))/(Options.canopyPoints-1);
        LADnorm(:,1) = 0;
    elseif Options.uniformLAD == 0
        zz = Canopy.z./Canopy.hN;
        LAD = exp(-((zz-0.4).^2) ./ (0.8));
        LAD(:,1) = 0;
        LAD = fliplr(LAD);
        LADnorm = [LAD(:,1)*0 LAD(:,1:end-1)]./sum(LAD,2);
    end
    Canopy.deltaLAI = Options.LAIRatio*LADnorm.*Weather.LAI;
    aLAI = [0*Canopy.deltaLAI(:,1) Canopy.deltaLAI];
    bLAI = [Canopy.deltaLAI 0*Canopy.deltaLAI(:,1)];
    dLAI = 0.5*(aLAI+bLAI); dLAI = dLAI(:,2:end);
    Canopy.LAI = fliplr(cumsum(fliplr(dLAI),2));
end

Canopy.LADProfile = Canopy.deltaLAI./Canopy.deltaZ; % Leaf area density [m-2 m-3]
for i= 1:1:length(Canopy.x)
    kB = @(psi) sqrt(Canopy.x(i).^2 + tan(psi).^2)./...
        (Canopy.x(i)+1.774.*(Canopy.x(i)+1.182).^(-0.733));
    tauB = @(psi) exp(-kB(psi).*Canopy.omega(i).*Weather.LAI(i));
    calltauD = @(psi) 2.*tauB(psi).*sin(psi).*cos(psi);
    tauD = integral(calltauD,0,pi/2);
    kD(i) = -log(tauD)./Canopy.omega(i)./Weather.LAI(i);
end
kD(isnan(kD)) = 0;
kD(isinf(kD)) = 0;
Canopy.kD = kD'; % Extinction coefficient for diffuse radiation [-]
Canopy.tauD = exp(-Canopy.kD.*Canopy.omega.*Canopy.deltaLAI); % Fraction of diffuse radiation 
% un-intercepted by the canopy layer [-] 

%% Initialising variables for the canopy microenvironment
Canopy.tAirProfileLag = zeros(height(Weather),(size(Canopy.deltaLAI,2))+1) + NaN;
Canopy.airQProfile = zeros(height(Weather),(size(Canopy.deltaLAI,2))+1) + NaN;
Canopy.airQProfileLag = zeros(height(Weather),(size(Canopy.deltaLAI,2))+1) + NaN;
Canopy.sunHeatSink = zeros(height(Weather),(size(Canopy.deltaLAI,2)));
Canopy.shadeHeatSink = zeros(height(Weather),(size(Canopy.deltaLAI,2)));
Canopy.sunTleaf = zeros(height(Weather),(size(Canopy.deltaLAI,2))) + Weather.tAir;
Canopy.shadeTleaf = zeros(height(Weather),(size(Canopy.deltaLAI,2))) + Weather.tAir;
Canopy.sunTleafLag = zeros(height(Weather),(size(Canopy.deltaLAI,2))) + NaN;
Canopy.shadeTleafLag = zeros(height(Weather),(size(Canopy.deltaLAI,2))) + NaN;
Canopy.tLeafProfile = zeros(height(Weather),(size(Canopy.deltaLAI,2))) + NaN;
Canopy.tLeafProfileLag  = zeros(height(Weather),(size(Canopy.deltaLAI,2))) + NaN;
Canopy.stemHeatSink = zeros(height(Weather),(size(Canopy.deltaLAI,2)+1));
Canopy.leafHeatSink = zeros(height(Weather),(size(Canopy.deltaLAI,2)+1));
end