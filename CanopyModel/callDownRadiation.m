%% Function to attenuate diffuse/scattered down radiation

function [Canopy] = callDownRadiation(EnergyOptions,Canopy,Options,layers)

%% Input
incident = sum(Canopy.down);
Error.threshold = 0.1;

if strcmp(Options.wavelength,'PAR')
    leafAbsorptivity = EnergyOptions.aPARLeaf; % Leaf absorptivity
    leafReflectivity = (1-leafAbsorptivity^0.5)/(1+leafAbsorptivity^0.5); % Leaf reflectivity
    leafTransmissivity = 1.0-leafAbsorptivity-leafReflectivity; % Leaf transmissivity
    soilAbsorptivity = EnergyOptions.aPARSoil;
    soilReflectivity = 1-EnergyOptions.aPARSoil;
    soilTransmissivity = 0;
elseif strcmp(Options.wavelength,'NIR')
    leafAbsorptivity = EnergyOptions.aNIRLeaf; % Leaf absorptivity Changed from 0.15
    leafReflectivity = (1-leafAbsorptivity^0.5)/(1+leafAbsorptivity^0.5); % Leaf reflectivity
    leafTransmissivity = 1.0-leafAbsorptivity-leafReflectivity; % Leaf transmissivity
    soilAbsorptivity = EnergyOptions.aNIRSoil;
    soilReflectivity = 1-EnergyOptions.aNIRSoil;
    soilTransmissivity = 0;
elseif strcmp(Options.wavelength,'long') % 
    leafAbsorptivity = EnergyOptions.epsilonLeaf; % Leaf absorptivity Changed from 0.96
    leafReflectivity = 0; %(1-leafAbsorptivity^0.5)/(1+leafAbsorptivity^0.5); % Leaf reflectivity Changed to 0 from Drewry
    leafTransmissivity = 1.0-leafAbsorptivity-leafReflectivity; % Leaf transmissivity
    soilAbsorptivity = EnergyOptions.epsilonSoil; % Changed from 0.8
    soilReflectivity = 1-EnergyOptions.epsilonSoil; % Changed from 0.2
    soilTransmissivity = 0;
end

Options.direction = 'down';

%% Initialize variables
Canopy.absorbed = zeros(1,length(Canopy.LAI));
Canopy.transmitted = zeros(1,length(Canopy.LAI));
Canopy.reflected = zeros(1,length(Canopy.LAI));

for layerLoop = 1:1:layers
    
    if strcmp(Options.orientation,'direct')
        k = sqrt(Canopy.x^2 + tand(Canopy.zenith).^2)/ ...
            (Canopy.x+1.774*(Canopy.x+1.182).^(-0.733));
    elseif strcmp(Options.orientation,'diffuse') % Canopy.zenith = Weather.zenith
        k = Canopy.kD;
    end
    
    extinction = exp(-k.*Canopy.omega.* ...
        [cumsum(Canopy.deltaLAI(2:end-layerLoop+1),'reverse') 0]); % Cumulative extinction coefficient at top of layer
    unintercepted = Canopy.down(1,end-layerLoop+1).*extinction; % Unintercepted radiation at the top of layer
    intercepted = [unintercepted(1) diff(unintercepted)]; % Cumulative intercepted radiation in each layer
    
    Canopy.absorbed(1,1:end-layerLoop+1) = Canopy.absorbed(1,1:end-layerLoop+1)+ ...
        [soilAbsorptivity*intercepted(1) leafAbsorptivity*intercepted(2:end)]; % Absorbed radiation
    Canopy.transmitted(1,1:end-layerLoop+1) = Canopy.transmitted(1,1:end-layerLoop+1)+ ...
        [soilTransmissivity*intercepted(1) leafTransmissivity*intercepted(2:end)]; % Transmitted radiation
    Canopy.reflected(1,1:end-layerLoop+1) = Canopy.reflected(1,1:end-layerLoop+1)+ ...
        [soilReflectivity*intercepted(1) leafReflectivity*intercepted(2:end)]; % Reflected radiation
end

%% Error analysis
Error.absolute = abs(incident-sum(Canopy.absorbed+Canopy.transmitted+Canopy.reflected));

if (Error.absolute >= Error.threshold)
    disp(strcat('Down radiation not converged at canopy ID = ',num2str(Canopy.ID), ...
        'for',Options.wavelength,'and for',Options.orientation))
end

end