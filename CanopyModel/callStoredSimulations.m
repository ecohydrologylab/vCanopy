if EnergyOptions.mixingSwitch == 0 && EnergyOptions.switch == 0
    savedcaProfile = Daily.savedcaProfile;
    savedwindProfile = Daily.savedwindProfile;
    savedeaProfile = Daily.savedeaProfile;
    savedtAirProfile = Daily.savedtAirProfile;
    savedRHProfile = Daily.savedRHProfile;
else
    savedwindProfile = zeros(length(Weather),size(Canopy(1).z,2)+1);
    savedcaProfile = zeros(length(Weather),size(Canopy(1).z,2)+1);
    savedeaProfile = zeros(length(Weather),size(Canopy(1).z,2)+1);
    savedtAirProfile = zeros(length(Weather),size(Canopy(1).z,2)+1);
    savedRHProfile = zeros(length(Weather),size(Canopy(1).z,2)+1);
end

if EnergyOptions.switch == 0
    savedtLeaf = Daily.savedtLeaf;
else
    savedtLeaf = zeros(length(Weather),2*size(Canopy(1).z,2));
end

if EnergyOptions.longwaveSwitch == 0 && EnergyOptions.switch == 0
    savedlongAbsorbed = Daily.savedlongAbsorbed;
    savedlongEmitted = Daily.savedlongEmitted;
else
    savedlongAbsorbed = zeros(length(Weather),size(Canopy(1).z,2));
    savedlongEmitted = zeros(length(Weather),size(Canopy(1).z,2));
end