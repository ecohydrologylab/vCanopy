function SSE = callOptimizer(X,Constants,Stomata,LeafMassFlux,LeafState,Photosynthesis,Weather)

    % Variable Vector   
    LeafMassFlux.aNet = X(1); LeafState.cbs = X(2); LeafState.ci = X(3);
    LeafState.gs = X(4); LeafState.cb = X(5);

    LeafState.g = LeafState.gs*LeafState.gb*Photosynthesis.Mbv/(LeafState.gs + LeafState.gb*Photosynthesis.Mbv); % Total leaf vapour conductance [mol m-2 sec-1]
    gv = LeafState.g;

    LeafState.eb = (LeafState.gs * LeafState.ei + LeafState.gb*Photosynthesis.Mbv * Weather.ea)...
        / (LeafState.gs + LeafState.gb*Photosynthesis.Mbv); % Vapour Pressure at leaf for boundary layer conductance [Pa]


    %% Preliminary Calculations for suitable Initial Calculations
    % Electron Transport rate
    PhiPS2 = 0.352 + 0.022 * LeafState.temperature - 3.4 * (LeafState.temperature)^2.0 / 10000.0;
    ThetaPS2 = Photosynthesis.theta;
    I = Weather.PAR*Constants.convert*PhiPS2/2;
    LeafMassFlux.J = (I + LeafMassFlux.jmax - ((I + LeafMassFlux.jmax)^(2.0) - 4.0 * ThetaPS2 * I * LeafMassFlux.jmax)^(0.5))/(2 * ThetaPS2);

    % Anet Calculation
    LeafMassFlux.vcLight = (LeafState.cbs*(1-Photosynthesis.x)*LeafMassFlux.J)/(3*LeafState.cbs+7*LeafState.GammaStar);
    LeafMassFlux.vcCO2 = (LeafState.cbs*LeafMassFlux.vcmax)/(LeafState.cbs+LeafState.Kc*(1+LeafState.obs/LeafState.Ko));
    LeafMassFlux.vc = min(LeafMassFlux.vcCO2,LeafMassFlux.vcLight);
    LeafMassFlux.aNet = (1-LeafState.GammaStar/LeafState.cbs)*LeafMassFlux.vc-LeafMassFlux.rd;

    % Oxygen Calculation
    LeafState.obs = Photosynthesis.alpha * LeafMassFlux.aNet /...
        (0.047 * Photosynthesis.gbs) + Weather.O2; % O2 conc. at bundle sheath [u moles mole-1]

    LeafMassFlux.vpLight = Photosynthesis.x*LeafMassFlux.J/2;
    LeafMassFlux.vpCO2 = min(LeafState.ci*LeafMassFlux.vpmax/(LeafState.ci+LeafState.Kp),100);
    LeafMassFlux.vp = min(LeafMassFlux.vpLight,LeafMassFlux.vpCO2);
    if LeafMassFlux.vcCO2 < LeafMassFlux.vcLight
        x1 = LeafMassFlux.vcmax;
        x2 = LeafState.Kc*(1+LeafState.obs/LeafState.Ko);            
    else
        x1 = (1-Photosynthesis.x)*LeafMassFlux.J/3;
        x2 = 7*LeafState.GammaStar/3;
    end
    
        %% Compensation point calculation
    LeafState.GammaStar = LeafState.gammaStar*LeafState.obs;
    LeafState.Gamma_C = (LeafState.GammaStar + LeafState.Kc*(1+LeafState.obs/LeafState.Ko)*LeafMassFlux.rd/LeafMassFlux.vcmax)/...
        (1-LeafMassFlux.rd/LeafMassFlux.vcmax); % Gamma for enzyme limited
    % end
    if LeafState.Gamma_C < 0
        LeafState.Gamma_C = 0;
    end

    GammaQuadratic = [(Photosynthesis.gbs) (Photosynthesis.gbs*(LeafState.Kp-LeafState.Gamma_C)+LeafMassFlux.vpmax-LeafMassFlux.rm)...
        -(LeafState.Gamma_C*LeafState.Kp*Photosynthesis.gbs+LeafMassFlux.rm*LeafState.Kp)];

    try
        Root = roots(GammaQuadratic);
        LeafState.Gamma = max(Root);
        if LeafState.Gamma < 0
            LeafState.Gamma = 2;
            disp('I am assuming the compensation point as 2')
        end
    catch
        LeafState.Gamma = 2;
    end

    % Function Vector
    SSE = ((LeafMassFlux.aNet - X(1))^2 + ...
    (LeafState.cbs - X(2))^2 + ... 
    (LeafState.ci - X(3))^2 + ...
    (LeafState.gs - X(4))^2 + ...
    (LeafState.cb - X(5))^2); 

        
end
