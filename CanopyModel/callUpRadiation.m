%% Function to attenuate diffuse/scattered down radiation

function [Canopy] = callUpRadiation(EnergyOptions,Canopy,Options,layers)

%% Input
incident = sum(Canopy.up(1:end-1));
skyLost = 0;

Error.threshold = 0.1;

if strcmp(Options.wavelength,'PAR')
    leafAbsorptivity = EnergyOptions.aPARLeaf; % Leaf absorptivity
    leafReflectivity = (1-leafAbsorptivity^0.5)/(1+leafAbsorptivity^0.5); % Leaf reflectivity
    leafTransmissivity = 1.0-leafAbsorptivity-leafReflectivity; % Leaf transmissivity
    soilAbsorptivity = EnergyOptions.aPARSoil;
    soilReflectivity = 1-EnergyOptions.aPARSoil;
    soilTransmissivity = 0;    
elseif strcmp(Options.wavelength,'NIR')
    leafAbsorptivity = EnergyOptions.aNIRLeaf; % Leaf absorptivity
    leafReflectivity = (1-leafAbsorptivity^0.5)/(1+leafAbsorptivity^0.5); % Leaf reflectivity
    leafTransmissivity = 1.0-leafAbsorptivity-leafReflectivity; % Leaf transmissivity
    soilAbsorptivity = EnergyOptions.aNIRSoil;
    soilReflectivity = 1-EnergyOptions.aNIRSoil;
    soilTransmissivity = 0;    
elseif strcmp(Options.wavelength,'long') % XXX check
    leafAbsorptivity = EnergyOptions.epsilonLeaf; % Leaf absorptivity
    leafReflectivity = 0;%(1-leafAbsorptivity^0.5)/(1+leafAbsorptivity^0.5); % Leaf reflectivity
    leafTransmissivity = 1.0-leafAbsorptivity-leafReflectivity; % Leaf transmissivity
    soilAbsorptivity = EnergyOptions.epsilonSoil;
    soilReflectivity = 1-EnergyOptions.epsilonSoil;
    soilTransmissivity = 0;    
end

Options.direction = 'up';

%% Initialize variables
Canopy.absorbed = zeros(1,length(Canopy.LAI));
Canopy.transmitted = zeros(1,length(Canopy.LAI));
Canopy.reflected = zeros(1,length(Canopy.LAI));

for layerLoop = 1:1:layers-1
    if(strcmp(Options.kDVariation,'Constant'))
        k = Canopy.kD;
    elseif (strcmp(Options.kDVariation,'Variable'))
        [k,tau] = callKDiffuse(Canopy.zenith,Canopy.x,Canopy.omega, ...
            [0 Canopy.deltaLAI(layerLoop+1:end) ],Options); % variable kD
    end
    
    extinction = exp(-k.*Canopy.omega.* ...
        cumsum([0 Canopy.deltaLAI(layerLoop+1:end) ])); % Cumulative extinction coefficient from canopy bottom
    unintercepted = Canopy.up(1,layerLoop).*extinction; % Unintercepted radiation at the top of layer
    intercepted = [0 -diff(unintercepted) ]; % Cumulative intercepted radiation in each layer
    
    if layerLoop==1
        Canopy.absorbed(1,layerLoop:end) = Canopy.absorbed(1,layerLoop:end)+ ...
            [soilAbsorptivity*intercepted(1) leafAbsorptivity*intercepted(layerLoop+1:end)]; % Absorbed radiation
        Canopy.transmitted(1,layerLoop:end) = Canopy.transmitted(1,layerLoop:end)+ ...
            [soilTransmissivity*intercepted(1) leafTransmissivity*intercepted(layerLoop+1:end)]; % Transmitted radiation
        Canopy.reflected(1,layerLoop:end) = Canopy.reflected(1,layerLoop:end)+ ...
            [soilReflectivity*intercepted(1) leafReflectivity*intercepted(layerLoop+1:end)]; % Reflected radiation
    else
        Canopy.absorbed(1,layerLoop:end) = Canopy.absorbed(1,layerLoop:end)+ ...
            [leafAbsorptivity*intercepted(1,1:end)]; % Absorbed radiation
        Canopy.transmitted(1,layerLoop:end) = Canopy.transmitted(1,layerLoop:end)+ ...
            [leafTransmissivity*intercepted(1,1:end)]; % Transmitted radiation
        Canopy.reflected(1,layerLoop:end) = Canopy.reflected(1,layerLoop:end)+ ...
            [leafReflectivity*intercepted(1,1:end)]; % Reflected radiation
    end
    
    skyLost = skyLost+unintercepted(end);
    Canopy.up(end) = Canopy.up(end)+unintercepted(end);
end

%% Error analysis
Error.absolute = abs(incident-sum(Canopy.absorbed+Canopy.transmitted+Canopy.reflected)- ...
    skyLost); % Adding sky lost

if (Error.absolute >= Error.threshold)
    disp(strcat('Up radiation not converged at canopy ID = ',num2str(Canopy.ID), ...
        'for',Options.wavelength,'and for',Options.orientation))
end

end