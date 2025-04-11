function [Photosynthesis,LeafState,LeafMassFlux,LeafEnergyFlux] = ...
    callOneTimeLeafInitialization(PhotosynthesisProperties,Weather,Initialize)
Photosynthesis = PhotosynthesisProperties; 
Photosynthesis = rmfield(Photosynthesis,{'kn','PcSeasonalVariation','linearDTD',...
    'Jmax_Vcmax','nitrogenDecay','decayStart','years'});
Photosynthesis = struct2table(Photosynthesis);
Photosynthesis = table2struct(Photosynthesis);

LeafState = table();
LeafState.flag = zeros(length(Photosynthesis),1);
LeafState.cbs = Initialize.cbs;
LeafState.obs = zeros(length(Photosynthesis),1)+Weather.O2;
LeafState.ci = Initialize.ci;
LeafState.cm = LeafState.ci*0.98;
LeafState.cb = Initialize.cb;
LeafState.gs = Initialize.gs;
LeafState.tLeaf = Initialize.tLeaf;
LeafState.ei = 1.2*(0.611*exp(17.502*LeafState.tLeaf./(240.97+LeafState.tLeaf)))*1000;
LeafState.eb = Initialize.eb;
LeafState.gbForced = zeros(length(Photosynthesis),1)+1;
LeafState.gbFree = zeros(length(Photosynthesis),1)+1;
LeafState.gb = zeros(length(Photosynthesis),1)+5;
LeafState.g = zeros(length(Photosynthesis),1)+0.4;
LeafState.gm = zeros(length(Photosynthesis),1)+0.4;
LeafState.Ko = zeros(length(Photosynthesis),1)+450;
LeafState.Kc = zeros(length(Photosynthesis),1)+80;
LeafState.Kp = zeros(length(Photosynthesis),1)+650;
LeafState.gammaStar = zeros(length(Photosynthesis),1)+0.000193;
LeafState.GammaStar = zeros(length(Photosynthesis),1)+0.143;
LeafState.Gamma = zeros(length(Photosynthesis),1)+1;
LeafState.Gamma_C = zeros(length(Photosynthesis),1)+20;
LeafState.Terror = zeros(length(Photosynthesis),1);
LeafState.Aerror = zeros(length(Photosynthesis),1);
LeafState.Cierror = zeros(length(Photosynthesis),1);
LeafState.Gserror = zeros(length(Photosynthesis),1);
LeafState.Cbserror = zeros(length(Photosynthesis),1);

LeafMassFlux = table();
LeafMassFlux.rd = zeros(length(Photosynthesis),1)+Photosynthesis(1).rd25.*Photosynthesis(1).vcmax25;
LeafMassFlux.rbs = LeafMassFlux.rd/2;
LeafMassFlux.rm = LeafMassFlux.rd/2;
LeafMassFlux.J = zeros(length(Photosynthesis),1)+Photosynthesis(1).jmax25;
LeafMassFlux.aNet = Initialize.aNet;
LeafMassFlux.tpu = zeros(length(Photosynthesis),1)+20;
LeafMassFlux.aGross = zeros(length(Photosynthesis),1)+2;
LeafMassFlux.ap = zeros(length(Photosynthesis),1);
LeafMassFlux.ac = zeros(length(Photosynthesis),1);
LeafMassFlux.aj = zeros(length(Photosynthesis),1);
LeafMassFlux.acCO2 = zeros(length(Photosynthesis),1);
LeafMassFlux.acLight = zeros(length(Photosynthesis),1);
LeafMassFlux.apCO2 = zeros(length(Photosynthesis),1);
LeafMassFlux.apLight = zeros(length(Photosynthesis),1);
LeafMassFlux.L_bs = zeros(length(Photosynthesis),1)+0.5;
LeafMassFlux.vp = zeros(length(Photosynthesis),1)+2;
LeafMassFlux.vc = zeros(length(Photosynthesis),1)+2;
LeafMassFlux.vpLight = zeros(length(Photosynthesis),1)+2;
LeafMassFlux.vpCO2 = zeros(length(Photosynthesis),1)+2;
LeafMassFlux.vcLight = zeros(length(Photosynthesis),1)+2;
LeafMassFlux.vcCO2 = zeros(length(Photosynthesis),1)+2;
LeafMassFlux.vpmax = zeros(length(Photosynthesis),1)+Photosynthesis(1).vpmax25;
LeafMassFlux.vcmax = zeros(length(Photosynthesis),1)+Photosynthesis(1).vcmax25;
LeafMassFlux.jmax = zeros(length(Photosynthesis),1)+Photosynthesis(1).jmax25;
LeafMassFlux.transpiration = zeros(length(Photosynthesis),1);

LeafEnergyFlux = table();
LeafEnergyFlux.H = zeros(length(Photosynthesis),1); % Sensible heat flux [W m-2]
LeafEnergyFlux.LE = zeros(length(Photosynthesis),1); % Latent heat flux [W m-2]
LeafEnergyFlux.storage = zeros(length(Photosynthesis),1); % Latent heat flux [W m-2]
LeafEnergyFlux.emission = zeros(length(Photosynthesis),1); % Long wave radiation flux emitted [W m-2]
LeafEnergyFlux.radiation = zeros(length(Photosynthesis),1); % Net radiation flux [W m-2]
LeafEnergyFlux.residual = 100*ones(length(Photosynthesis),1); % Error radiation flux [W m-2]

LeafState = table2struct(LeafState);
LeafMassFlux = table2struct(LeafMassFlux);
LeafEnergyFlux = table2struct(LeafEnergyFlux);



end