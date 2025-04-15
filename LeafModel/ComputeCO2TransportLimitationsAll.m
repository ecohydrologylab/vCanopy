% Total leaf conductance for water vapour
LeafState.g = LeafState.gs.*LeafState.gb.*Stomata.Mb./(LeafState.gs + LeafState.gb.*Stomata.Mb); % Total leaf vapour conductance [mol m-2 sec-1]
gv = LeafState.g; % Total leaf vapour conductance [mol m-2 sec-1]

% Boudnary layer water vapour pressure
LeafState.eb = (LeafState.gs ./ convert .* LeafState.ei + LeafState.gb.*Stomata.Mb  ./ convert .* Weather.ea)...
    ./ (LeafState.gs ./ convert + LeafState.gb.*Stomata.Mb ./ convert); % Vapour Pressure at leaf for boundary layer conductance [Pa]

%% Photosynthesis module: Preliminary Calculations for suitable Initial Calculations
% Electron Transport rate
PhiPS2 = 0.352 + 0.022 * LeafState.tLeaf - 3.4 .* (LeafState.tLeaf).^2.0 / 10000.0;
ThetaPS2 = Photosynthesis.theta;
I = Weather.PAR.*Constants.convert.*PhiPS2/2; % Radiation absorbed on leaf surface (PPFD) [u mol m-2 sec-1]
LeafMassFlux.J = (I + LeafMassFlux.jmax - ((I + LeafMassFlux.jmax).^(2.0) - 4.0 .* ThetaPS2 .* I .* LeafMassFlux.jmax).^(0.5))...
    ./(2 .* ThetaPS2); % Total chain electron transport rate [u mol m-2 sec-1]

% Rubisco carboxylation
LeafMassFlux.vcLight = (LeafState.cbs.*(1-Photosynthesis.x).*LeafMassFlux.J)...
    ./(3*LeafState.cbs+7*LeafState.GammaStar); % Light limited rubisco carboxylation [mu mol m-2 sec-1]
LeafMassFlux.vcCO2 = (LeafState.cbs.*LeafMassFlux.vcmax)....
    /(LeafState.cbs+LeafState.Kc.*(1+LeafState.obs./LeafState.Ko)); % CO2 limited rubisco carboxylation [u mol m-2 sec-1]
LeafMassFlux.vc = min(LeafMassFlux.vcCO2,LeafMassFlux.vcLight); % Rubisco carboxylation rate [u mol m-2 sec-1]

% Oxygen Calculation
LeafState.obs = Photosynthesis.alpha .* LeafMassFlux.aNet ./...
    (0.047 .* Photosynthesis.gbs.* 1000) + Weather.O2; % O2 conc. at bundle sheath [u moles mole-1]
LeafMassFlux.vpLight = Photosynthesis.x.*LeafMassFlux.J/2; % Light limited PEP carboxylation [mu mol m-2 sec-1]
LeafMassFlux.vpCO2 = min(LeafState.ci.*LeafMassFlux.vpmax...
    ./(LeafState.ci+LeafState.Kp),100); % CO2 limited PEP carboxylation [mu mol m-2 sec-1]
LeafMassFlux.vp = min(LeafMassFlux.vpLight,LeafMassFlux.vpCO2); % PEP carboxylation rate [mu mol m-2 sec-1]

% Deciding Anet limiting factors
if LeafMassFlux.vcCO2 < LeafMassFlux.vcLight
    a1 = LeafMassFlux.vcmax;
    b1 = LeafState.Kc.*(1+LeafState.obs./LeafState.Ko);            
else
    a1 = (1-Photosynthesis.x).*LeafMassFlux.J/3;
    b1 = 7*LeafState.GammaStar/3;
end

% Sub-limitation calculations for photosynthesis
LeafMassFlux.acCO2 = (1.0 - LeafState.gammaStar .* LeafState.obs ./ LeafState.cbs).*LeafMassFlux.vcCO2-...
    LeafMassFlux.rd; % CO2 limited rubisco Photosynthesis rate [mu mol m-2 sec-1]
LeafMassFlux.acLight = (1.0 - LeafState.gammaStar .* LeafState.obs ./ LeafState.cbs).*LeafMassFlux.vcLight-...
    LeafMassFlux.rd; % Light limited rubisco Photosynthesis rate [mu mol m-2 sec-1]
LeafMassFlux.apCO2 = (1.0 - LeafState.gammaStar .* LeafState.obs ./ LeafState.cbs).*LeafMassFlux.vpCO2-...
    LeafMassFlux.rd; % CO2 limited PEP Photosynthesis rate [mu mol m-2 sec-1]
LeafMassFlux.apLight = (1.0 - LeafState.gammaStar .* LeafState.obs ./ LeafState.cbs).*LeafMassFlux.vpLight-...
    LeafMassFlux.rd; % Light limited PEP Photosynthesis rate [mu mol m-2 sec-1]

% Compensation point calculation
LeafState.GammaStar = LeafState.gammaStar.*LeafState.obs; % Bundle sheat [CO2] at Anet  = 0 [ppm]
LeafState.Gamma_C = (LeafState.GammaStar + LeafState.Kc.*(1+LeafState.obs./LeafState.Ko).*LeafMassFlux.rd./LeafMassFlux.vcmax)./...
    (1-LeafMassFlux.rd./LeafMassFlux.vcmax); % CO2 compensation point without considering dark respiration [ppm]
if LeafState.Gamma_C < 0
    LeafState.Gamma_C = 0;
end

GammaQuadratic = [(Photosynthesis.gbs) (Photosynthesis.gbs.*(LeafState.Kp-LeafState.Gamma_C)+LeafMassFlux.vpmax-LeafMassFlux.rm)...
    -(LeafState.Gamma_C.*LeafState.Kp.*Photosynthesis.gbs+LeafMassFlux.rm.*LeafState.Kp)];

Gamma = zeros(length(Weather.PAR),1);
for i=1:1:size(GammaQuadratic,1)
    try
        Root = roots(GammaQuadratic(i,:));
        Gamma(i) = max(Root);
        if Gamma(i) < 0
            Gamma(i) = 2;
            disp('I am assuming the compensation point as 2')
        end
    catch
        Gamma(i) = 2;
    end
end 
LeafState.Gamma = Gamma; % CO2 compensation point after considering dark respirationn [ppm]