function [Canopy] = callInitializeCanopy(Weather,Canopy)

Canopy.PARDirectAbsorbed = zeros(height(Weather),size(Canopy.z,2));
Canopy.PARDirectTransmitted = zeros(height(Weather),size(Canopy.z,2));
Canopy.PARDirectReflected = zeros(height(Weather),size(Canopy.z,2));
Canopy.PARDiffuseAbsorbed = zeros(height(Weather),size(Canopy.z,2));
Canopy.PARDiffuseTransmitted = zeros(height(Weather),size(Canopy.z,2));
Canopy.PARDiffuseReflected = zeros(height(Weather),size(Canopy.z,2));
Canopy.PARScatteredAbsorbed = zeros(height(Weather),size(Canopy.z,2));
Canopy.PARScatteredUp = zeros(height(Weather),size(Canopy.z,2));
Canopy.sunPARAbsorbed = zeros(height(Weather),size(Canopy.z,2));
Canopy.shadePARAbsorbed = zeros(height(Weather),size(Canopy.z,2));

Canopy.NIRDirectAbsorbed = zeros(height(Weather),size(Canopy.z,2));
Canopy.NIRDirectTransmitted = zeros(height(Weather),size(Canopy.z,2));
Canopy.NIRDirectReflected = zeros(height(Weather),size(Canopy.z,2));
Canopy.NIRDiffuseAbsorbed = zeros(height(Weather),size(Canopy.z,2));
Canopy.NIRDiffuseTransmitted = zeros(height(Weather),size(Canopy.z,2));
Canopy.NIRDiffuseReflected = zeros(height(Weather),size(Canopy.z,2));
Canopy.NIRScatteredAbsorbed = zeros(height(Weather),size(Canopy.z,2));
Canopy.NIRScatteredUp = zeros(height(Weather),size(Canopy.z,2));
Canopy.sunNIRAbsorbed = zeros(height(Weather),size(Canopy.z,2));
Canopy.shadeNIRAbsorbed = zeros(height(Weather),size(Canopy.z,2));

Canopy.longEmitted  = zeros(height(Weather),size(Canopy.z,2));
Canopy.longAbsorbed = zeros(height(Weather),size(Canopy.z,2));
Canopy.longUp = zeros(height(Weather),size(Canopy.z,2));
Canopy.sunLongAbsorbed = zeros(height(Weather),size(Canopy.z,2));
Canopy.shadeLongAbsorbed = zeros(height(Weather),size(Canopy.z,2));

Canopy.sunFraction = zeros(height(Weather),size(Canopy.z,2));
Canopy.shadeFraction = zeros(height(Weather),size(Canopy.z,2));

Canopy.sunTemperature = repmat(Weather.temperature,[1,size(Canopy.z,2)]);
Canopy.shadeTemperature = repmat(Weather.temperature,[1,size(Canopy.z,2)]);

Canopy.aNet = zeros(height(Weather),length(Canopy.z(1,:)));
Canopy.latentHeat = zeros(height(Weather),length(Canopy.z(1,:)));
Canopy.sensibleHeat = zeros(height(Weather),length(Canopy.z(1,:)));
Canopy.biomassHeat = zeros(height(Weather),length(Canopy.z(1,:)));
Canopy.taud = zeros(height(Weather),length(Canopy.z(1,:)));

Canopy.FHs = zeros(height(Weather),1);
Canopy.FLE = zeros(height(Weather),1);
Canopy.FCO2 = zeros(height(Weather),1);

Canopy.windProfile = zeros(height(Weather),size(Canopy.deltaLAI,2)+1);
Canopy.temperatureProfile = zeros(height(Weather),(size(Canopy.deltaLAI,2))+1);
Canopy.airtempLag1 = zeros(height(Weather),(size(Canopy.deltaLAI,2))+1)+NaN;
Canopy.eaProfile = zeros(height(Weather),(size(Canopy.deltaLAI,2))+1);
Canopy.RHProfile = zeros(height(Weather),(size(Canopy.deltaLAI,2))+1);
Canopy.caProfile = zeros(height(Weather),(size(Canopy.deltaLAI,2))+1);

Canopy.flag = ones(height(Weather),2*(size(Canopy.deltaLAI,2)));
Canopy.Transpiration = zeros(height(Weather),2*(size(Canopy.deltaLAI,2)));
Canopy.LW_leaf = zeros(height(Weather),2*(size(Canopy.deltaLAI,2)));
Canopy.ANET = zeros(height(Weather),2*(size(Canopy.deltaLAI,2)));
Canopy.gs = zeros(height(Weather),2*(size(Canopy.deltaLAI,2)));
Canopy.LatentHeat = zeros(height(Weather),2*(size(Canopy.deltaLAI,2)));
Canopy.SensibleHeat = zeros(height(Weather),2*(size(Canopy.deltaLAI,2)));
Canopy.Residual = zeros(height(Weather),2*(size(Canopy.deltaLAI,2)));
Canopy.LAIflag = zeros(height(Weather),2*(size(Canopy.deltaLAI,2)));


end