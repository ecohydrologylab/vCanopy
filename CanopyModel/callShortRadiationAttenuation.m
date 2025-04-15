function [Soil,Canopy] = callShortRadiationAttenuation(EnergyOptions,Weather,Soil,Canopy)

%% Compute direct PAR radiation attenuation
Options.isotropy = EnergyOptions.isotropy;
Options.wavelength = 'PAR';
Options.orientation = 'direct';
Canopy.down = zeros(1,size(Canopy.z,2));
Canopy.down(end) = Weather.PARDirect;
Canopy.zenith = Weather.zenith;
layers = 1;

[Canopy] = callDownRadiation(EnergyOptions,Canopy,Options,layers);

Canopy.PARDirectAbsorbed = Canopy.absorbed;
Canopy.PARDirectTransmitted = Canopy.transmitted;
Canopy.PARDirectReflected = Canopy.reflected;

%% Compute diffuse PAR radiation attenuation

Options.wavelength = 'PAR';
Options.orientation = 'diffuse';
Canopy.down = zeros(1,size(Canopy.z,2));
Canopy.down(end) = Weather.PARDiffuse;
layers = 1;

[Canopy] = callDownRadiation(EnergyOptions,Canopy,Options,layers);

Canopy.PARDiffuseAbsorbed = Canopy.absorbed;
Canopy.PARDiffuseTransmitted = Canopy.transmitted;
Canopy.PARDiffuseReflected = Canopy.reflected;

%% Compute scattered PAR radiation attenuation

Canopy.scatteredDown = [Canopy.PARDirectTransmitted(2:end)+ ...
    Canopy.PARDiffuseTransmitted(2:end) 0]; % Shift as transmitted is bottom of layer
Canopy.scatteredUp = Canopy.PARDirectReflected+Canopy.PARDiffuseReflected;
Canopy.scatteredAbsorbed = zeros(1,size(Canopy.z,2));
Options.wavelength = 'PAR';
Options.orientation = 'diffuse';

[Canopy] = callScatteredRadiation(EnergyOptions,Canopy,Options);

Canopy.PARScatteredAbsorbed = Canopy.scatteredAbsorbed;
Canopy.PARScatteredUp = Canopy.scatteredUp;

%% Compute direct NIR radiation attenuation

Options.wavelength = 'NIR';
Options.orientation = 'direct';
Canopy.down = zeros(1,size(Canopy.z,2));
Canopy.down(end) = Weather.NIRDirect;
layers = 1;

[Canopy] = callDownRadiation(EnergyOptions,Canopy,Options,layers);

Canopy.NIRDirectAbsorbed = Canopy.absorbed;
Canopy.NIRDirectTransmitted = Canopy.transmitted;
Canopy.NIRDirectReflected = Canopy.reflected;

%% Compute diffuse NIR radiation attenuation

Options.wavelength = 'NIR';
Options.orientation = 'diffuse';
Canopy.down = zeros(1,size(Canopy.z,2));
Canopy.down(end) = Weather.NIRDiffuse;
layers = 1;

[Canopy] = callDownRadiation(EnergyOptions,Canopy,Options,layers);
Canopy.NIRDiffuseAbsorbed = Canopy.absorbed;
Canopy.NIRDiffuseTransmitted = Canopy.transmitted;
Canopy.NIRDiffuseReflected = Canopy.reflected;

%% Compute scattered NIR radiation attenuation

Canopy.scatteredDown = [Canopy.NIRDirectTransmitted(2:end)+ ...
    Canopy.NIRDiffuseTransmitted(2:end) 0]; % Shift as transmitted is bottom of layer
Canopy.scatteredUp = Canopy.NIRDirectReflected+Canopy.NIRDiffuseReflected;
Canopy.scatteredAbsorbed = zeros(1,size(Canopy.z,2));
Options.wavelength = 'NIR';
Options.orientation = 'diffuse';

[Canopy] = callScatteredRadiation(EnergyOptions,Canopy,Options);

Canopy.NIRScatteredAbsorbed = Canopy.scatteredAbsorbed;
Canopy.NIRScatteredUp = Canopy.scatteredUp;

%% Error
Error.PARAbsolute = Weather.PARDirect+Weather.PARDiffuse- ...
    sum(Canopy.PARDirectAbsorbed+Canopy.PARDiffuseAbsorbed+ ...
    Canopy.PARScatteredAbsorbed)-Canopy.PARScatteredUp(end); % PAR radiation
Error.NIRAbsolute = Weather.NIRDirect+Weather.NIRDiffuse- ...
    sum(Canopy.NIRDirectAbsorbed+Canopy.NIRDiffuseAbsorbed+ ...
    Canopy.NIRScatteredAbsorbed)-Canopy.NIRScatteredUp(end); % NIR radiation

%% Compute air and soil radiation

% Air.PAROut = Canopy.PARScatteredUp(end);
Soil.PARDirectAbsorbed = Canopy.PARDirectAbsorbed(1);
Soil.PARDiffuseAbsorbed = Canopy.PARDiffuseAbsorbed(1);
Soil.PARScatteredAbsorbed = Canopy.PARScatteredAbsorbed(1);
Soil.PARAbsorbed = Canopy.PARDirectAbsorbed(1)+Canopy.PARDiffuseAbsorbed(1)+ ...
    Canopy.PARScatteredAbsorbed(1);

% Air.NIROut = Canopy.NIRScatteredUp(end);
Soil.NIRDirectAbsorbed = Canopy.NIRDirectAbsorbed(1);
Soil.NIRDiffuseAbsorbed = Canopy.NIRDiffuseAbsorbed(1);
Soil.NIRScatteredAbsorbed = Canopy.NIRScatteredAbsorbed(1);
Soil.NIRAbsorbed = Canopy.NIRDirectAbsorbed(1)+Canopy.NIRDiffuseAbsorbed(1)+ ...
    Canopy.NIRScatteredAbsorbed(1);

%% PAR absorbed sun-shade model

Canopy.sunPARAbsorbed = Canopy.PARDirectAbsorbed+Canopy.sunFraction.* ...
    (Canopy.PARDiffuseAbsorbed+Canopy.PARScatteredAbsorbed);
Canopy.shadePARAbsorbed = Canopy.shadeFraction.* ...
    (Canopy.PARDiffuseAbsorbed+Canopy.PARScatteredAbsorbed);

%% NIR absorbed sun-shade model

Canopy.sunNIRAbsorbed = Canopy.NIRDirectAbsorbed+Canopy.sunFraction.* ...
    (Canopy.NIRDiffuseAbsorbed+Canopy.NIRScatteredAbsorbed);
Canopy.shadeNIRAbsorbed = Canopy.shadeFraction.* ...
    (Canopy.NIRDiffuseAbsorbed+Canopy.NIRScatteredAbsorbed);

%% Remove duplicate fields in structure

Canopy = rmfield(Canopy,{'down','up','scatteredUp','scatteredDown', ...
    'scatteredAbsorbed','absorbed','transmitted','reflected', ...
    'zenith','PARDirectAbsorbed','PARDiffuseAbsorbed','PARScatteredAbsorbed',...
    'NIRDirectAbsorbed','NIRDiffuseAbsorbed','NIRScatteredAbsorbed',...
    'PARDirectTransmitted','PARDirectReflected','PARDiffuseTransmitted',...
    'PARDiffuseReflected','NIRDirectTransmitted','NIRDirectReflected',...
    'NIRDiffuseTransmitted','NIRDiffuseReflected'});

end