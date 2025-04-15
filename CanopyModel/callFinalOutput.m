%% RequiredOutput
CanopyOut.sunFlag = [LeafState(1:nLayer).flag];
CanopyOut.sunFraction = Canopy.sunFraction*nanFlag;
CanopyOut.sunPARAbsorbed = Canopy.sunPARAbsorbed*nanFlag;
CanopyOut.sunNIRAbsorbed = Canopy.sunNIRAbsorbed*nanFlag;
CanopyOut.sunLongAbsorbed = Canopy.sunLongAbsorbed*nanFlag;
CanopyOut.sunGs = [LeafState(1:nLayer).gs];
CanopyOut.sunCi = [LeafState(1:nLayer).ci];
CanopyOut.sunCb = [LeafState(1:nLayer).cb];
CanopyOut.sunTleaf = Canopy.sunTleaf;
CanopyOut.sunCbs = [LeafState(1:nLayer).cbs];
CanopyOut.sunAnet = [LeafMassFlux(1:nLayer).aNet];
CanopyOut.sunTranspiration = [LeafMassFlux(1:nLayer).transpiration];
CanopyOut.sunLE = [LeafEnergyFlux(1:nLayer).LE];
CanopyOut.sunH = [LeafEnergyFlux(1:nLayer).H];
CanopyOut.sunResidual = [LeafEnergyFlux(1:nLayer).residual];
CanopyOut.sunEb = [LeafState(1:nLayer).eb];
CanopyOut.sunGammaStar = [LeafState(1:nLayer).GammaStar];

CanopyOut.shadeFlag = [LeafState(nLayer+1:end).flag];
CanopyOut.shadeFraction = Canopy.shadeFraction*nanFlag;
CanopyOut.shadePARAbsorbed = Canopy.shadePARAbsorbed*nanFlag;
CanopyOut.shadeNIRAbsorbed = Canopy.shadeNIRAbsorbed*nanFlag;
CanopyOut.shadeLongAbsorbed = Canopy.shadeLongAbsorbed*nanFlag;
CanopyOut.shadeGs = [LeafState(nLayer+1:end).gs];
CanopyOut.shadeCi = [LeafState(nLayer+1:end).ci];
CanopyOut.shadeCb = [LeafState(nLayer+1:end).cb];
CanopyOut.shadeTleaf = Canopy.shadeTleaf;
CanopyOut.shadeCbs = [LeafState(nLayer+1:end).cbs];
CanopyOut.shadeAnet = [LeafMassFlux(nLayer+1:end).aNet];
CanopyOut.shadeTranspiration = [LeafMassFlux(nLayer+1:end).transpiration];
CanopyOut.shadeLE = [LeafEnergyFlux(nLayer+1:end).LE];
CanopyOut.shadeH = [LeafEnergyFlux(nLayer+1:end).H];
CanopyOut.shadeResidual = [LeafEnergyFlux(nLayer+1:end).residual];
CanopyOut.shadeEb = [LeafState(nLayer+1:end).eb];
CanopyOut.shadeGammaStar = [LeafState(nLayer+1:end).GammaStar];

CanopyOut.tLeafProfile = Canopy.tLeafProfile*nanFlag;
CanopyOut.airQProfile = Canopy.airQProfile*nanFlag;
CanopyOut.aNet = Canopy.aNet*nanFlag;
CanopyOut.LE = Canopy.LE*nanFlag;
CanopyOut.H = Canopy.H*nanFlag;
CanopyOut.longEmitted = Canopy.longEmitted*nanFlag;
CanopyOut.longAbsorbed = Canopy.longAbsorbed*nanFlag;
CanopyOut.longUp = Canopy.longUp*nanFlag;
CanopyOut.PARScatteredUp = Canopy.PARScatteredUp*nanFlag;
CanopyOut.NIRScatteredUp = Canopy.NIRScatteredUp*nanFlag;
CanopyOut.stemHeatSink = Canopy.stemHeatSink*nanFlag;
CanopyOut.leafHeatSink = Canopy.leafHeatSink*nanFlag;
CanopyOut.windProfile = Canopy.windProfile*nanFlag;
CanopyOut.caProfile = Canopy.caProfile*nanFlag;
CanopyOut.tAirProfile = Canopy.tAirProfile*nanFlag;
CanopyOut.eaProfile = Canopy.eaProfile*nanFlag;
CanopyOut.RHProfile = Canopy.RHProfile*nanFlag;
CanopyOut.FCO2 = Canopy.FCO2*nanFlag;
CanopyOut.FLE = Canopy.FLE*nanFlag;
CanopyOut.FHs = Canopy.FHs*nanFlag;
CanopyOut.soiltemperature = Soil.temperature*nanFlag;

logFile = "";

try
    if Error.loop >= Error.maxLoop && ~isnan(nanFlag)
        disp(strcat("Aborting canopy microenvironment for weather, ",num2str(Weather.hour)," hr; Input row = ",num2str(loop)))
        logFile = strcat(logFile,"Microenvironment_convergence_loop,",num2str(loop)," \n");
    end
    
    if sum(isnan(CanopyOut.caProfile(2:end))) || sum(isnan(CanopyOut.tAirProfile(2:end))) ...
            || sum(isnan(CanopyOut.eaProfile(2:end)))
        disp(strcat("Aborting canopy microenvironment for weather, ",num2str(Weather.hour)," hr; Input row = ",num2str(loop)))
        logFile = strcat(logFile,"Microenvironment_complex_loop,",num2str(loop)," \n");
    end
    
    for i = 1:1:length(Canopy.z)
        if CanopyOut.sunFlag(i) == 1
            disp(strcat("Leaf did not converge at ,",num2str(Weather.hour),", hr loop ,",num2str(loop),", for leafLoop ,",num2str(i)));
            logFile = strcat(logFile,"Leaf not converge at weather loop,",num2str(loop),", leafLoop,",num2str(i)," \n");
        end
        if CanopyOut.shadeFlag(i) == 1
            disp(strcat("Leaf did not converge at ,",num2str(Weather.hour),", hr loop ,",num2str(loop),", for leafLoop ,",num2str(i)));
            logFile = strcat(logFile,"Leaf not converge at weather loop,",num2str(loop),", leafLoop,",num2str(i+length(Canopy.z))," \n");
        end
        
    end
catch
    logFile = strcat(logFile,"\n Problem at ",num2str(loop)," \n");
end