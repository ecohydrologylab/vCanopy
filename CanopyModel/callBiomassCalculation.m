
%% Calculation of the developmental stage
DVS = ones(height(Weather),1);
index = Weather.GDD/Options.GDDpre < 1;
DVS(index) = 0.5*(Weather.GDD(index)/Options.GDDpre);
index = Weather.GDD/Options.GDDpre >= 1;
DVS(index) = 0.5 + 0.5*(Weather.GDD(index) - Options.GDDpre)/Options.GDDpost;

%% Calculation of slecific leaf area
Canopy.SLA = 13.37 - 5.87*log(DVS); % [m2(leaf) kg-1] Meter square per kg of dry weight

%% Canopy Dry Leaf Mass and Dry Stem Mass
Canopy.dryLM = Canopy.deltaLAI./Canopy.SLA; % [kg m-2(ground)] leaf dry Weight per unit ground area
Canopy.drySM = 6.05*Canopy.dryLM.^0.827; % [kg m-2(fround)] stem dry weight per unit ground area

