%% Function to create canopy input
function [Soil,Canopy] = callSoilCanopyAirInput(Weather,Options)

%% Soil properties
Soil = table;
Soil.temperature = Weather.tAir+1;
Soil.wind = zeros(height(Weather),1); % Wind speed at soil [m s-1]
Soil.RH = 85+zeros(height(Weather),1); % Relative humidity [%]
Soil.roughnessHeight = 0.005*ones(height(Weather),1);
Soil.CO2Flux = 1*ones(height(Weather),1);
Soil.LE = 1*ones(height(Weather),1);
Soil.H = 1*ones(height(Weather),1);
Soil.G = 1*ones(height(Weather),1);
Soil.netRadiation = zeros(height(Weather),1);

%% Canopy properties
Canopy = table;
if Options.plant == 3
    Canopy.hN = 1*ones(height(Weather),1); % End height [m]
elseif Options.plant == 4
    try
        Canopy.hN = Weather.Height;
    catch
        Canopy.hN = 2.5*ones(height(Weather),1); % End height [m]
    end
end
Canopy.ID = [1:1:height(Weather)]';
deltaZ = Canopy.hN/(Options.canopyPoints-1);
Canopy.deltaZ = [-1*deltaZ/2,repmat(deltaZ,[1,Options.canopyPoints-1])];
Canopy.z = cumsum(Canopy.deltaZ,2);
Canopy.deltaZ(:,1) = NaN;

Canopy.lm = 0.1367 * Canopy.hN;
Canopy.cd = 0.2*ones(height(Weather),1);
Canopy.x = Options.X*ones(height(Weather),1); % Leaf angle distribution factor Soybean
Canopy.omega = Options.omega*ones(height(Weather),1); % Clumping factor

if Options.canopyPoints == 2 % For big leaf model
    Canopy.deltaLAI = [Weather.LAI*0 Options.LAIRatio*Weather.LAI]; % Add soil layer
    Canopy.LAI = Canopy.deltaLAI;
else
    
    if Options.uniformLAD == 1
        LADnorm = ones(size(Canopy.z,1),size(Canopy.z,2))/(Options.canopyPoints-1);
        LADnorm(:,1) = 0;
    elseif Options.uniformLAD == 0
        zz = Canopy.z./Canopy.hN;
        if Options.plant == 4
            LAD = exp(-((zz-0.4).^2) ./ (0.8));
            LAD(:,1) = 0;
            LAD = fliplr(LAD);
            LADnorm = [LAD(:,1)*0 LAD(:,1:end-1)]./sum(LAD,2);
        elseif Options.plant == 3
            % Observed data for 10 layers
            Obs_z_norm = [5:10:95]/100;
            Obs_LAI_norm = [0.0000, 0.0993, 0.0880, 0.0817, 0.1005, 0.0894,...
                0.1094, 0.1532, 0.2162, 0.0624];
            % Interpolation to the new grid
            LADnorm = interp1(Obs_z_norm, Obs_LAI_norm, zz, 'linear', 'extrap');
            LADnorm(LADnorm<0) = 0;
            LADnorm = LADnorm./sum(LADnorm,2);
        end
    end
    Canopy.deltaLAI = Options.LAIRatio*LADnorm.*Weather.LAI;
    aLAI = [0*Canopy.deltaLAI(:,1) Canopy.deltaLAI];
    bLAI = [Canopy.deltaLAI 0*Canopy.deltaLAI(:,1)];
    dLAI = 0.5*(aLAI+bLAI); dLAI = dLAI(:,2:end);
    Canopy.LAI = fliplr(cumsum(fliplr(dLAI),2));
end

Canopy.LADProfile = Canopy.deltaLAI./Canopy.deltaZ;
if strcmp(Options.kDVariation,"Constant")
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
    Canopy.kD = kD';
    Canopy.tauD = exp(-Canopy.kD.*Canopy.omega.*Canopy.deltaLAI);
else
    for i= 1:1:length(Canopy.x)
        for j = 2:1:size(Canopy.deltaLAI,2)
            kB = @(psi) sqrt(Canopy.x(i).^2 + tan(psi).^2)./...
                (Canopy.x(i)+1.774.*(Canopy.x(i)+1.182).^(-0.733));
            tauB = @(psi) exp(-kB(psi).*Canopy.omega(i).*Canopy.deltaLAI(i,j));
            calltauD = @(psi) 2.*tauB(psi).*sin(psi).*cos(psi);
            tauD = integral(calltauD,0,pi/2);
            kD(i,j) = -log(tauD)./Canopy.omega(i)./Canopy.deltaLAI(i,j);
        end
    end
    kD(isnan(kD)) = 0;
    kD(isinf(kD)) = 0;
    Canopy.kD = kD;
    Canopy.tauD = exp(-Canopy.kD.*Canopy.omega.*Canopy.deltaLAI); 
    Canopy.tauD(isnan(Canopy.tauD)) = 1;
end

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