function [LeafState] = callStomatalConductance(Stomata,LeafState,LeafMassFlux,Weather)

%% Ball Berry equation for stomatal conductance [mol m-2 sec-1]
LeafState.gs = max(Stomata.intercept, Stomata.intercept+...
    Stomata.fsv*Stomata.slope*LeafMassFlux.aNet*(LeafState.eb/LeafState.ei)/...
    (LeafState.cb-LeafState.Gamma)); % Stomatal conductance for vapour
%     [mol m-2 s-1]  /Weather.ca)

%% Mass balance equation for ci [ppm]
LeafState.ci = LeafState.cb - 1.6*LeafMassFlux.aNet/LeafState.gs; % Intercellular concentration of CO2 in air corrected for solubility relative to 25 degree Celcius [ppm]

end