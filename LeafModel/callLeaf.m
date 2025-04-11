%% Function to compute leaf photosynthesis
function [LeafMassFlux,LeafEnergyFlux,LeafState] = callLeaf(Constants,...
    LeafBoundaryLayer,EnergyOptions,Photosynthesis,Stomata,Weather,CanopyLayer,...
    LeafState,LeafMassFlux,LeafEnergyFlux,loop)

%% Convergence criteria
Error.ciMax = 0.1; % [ppm]
Error.gsMax = 0.01; % [mol m-2 s-1]
Error.aNetMax = 0.05; % [u mol m-2 s-1]
Error.tLeafMax = 0.05; % [Celcius]
Error.loop = 1;
Error.loopMax = 500;
Error.relax = 0.5; % Relaxation factor
Error.gs = 100;
Error.aNet = 100;
Error.ci = 100;
Error.tLeaf = 100;

PreviousLeafMassFlux = LeafMassFlux;
PreviousLeafState = LeafState;

%% Gauss Seidal Method for leaf solution with successive relaxation
while ((Error.loop < Error.loopMax) && (Error.ci >= Error.ciMax || ...
        Error.gs >= Error.gsMax || Error.aNet >= Error.aNetMax...
        || Error.tLeaf >= Error.tLeafMax))
    
    % Compute photosynthesis
    [LeafMassFlux,LeafState] = callC4Photosynthesis(Constants,Photosynthesis,Weather,LeafState,LeafMassFlux,loop);

    LeafMassFlux.aNet =  Error.relax*LeafMassFlux.aNet + (1-Error.relax)*PreviousLeafMassFlux.aNet; % Apply relaxation to prevent oscillation
    LeafState.cbs = Error.relax*LeafState.cbs + (1-Error.relax)*PreviousLeafState.cbs; % Apply relaxation to prevent oscillation
    
    % Compute Boundary Layer Conductance
    [LeafState] = callBoundaryLayerConductance(LeafBoundaryLayer,Stomata,Weather,LeafState,LeafMassFlux);
    LeafState.cb = Error.relax*LeafState.cb + (1-Error.relax)*PreviousLeafState.cb; % Apply relaxation to prevent oscillation
    
    % Compute stomatal conductance
    [LeafState] = callStomatalConductance(Stomata,LeafState,LeafMassFlux,Weather);
    LeafState.ci = Error.relax*LeafState.ci + (1-Error.relax)*PreviousLeafState.ci; % Apply relaxation to prevent oscillation
    LeafState.gs = Error.relax*LeafState.gs + (1-Error.relax)*PreviousLeafState.gs; % Apply relaxation to prevent oscillation
    
    % Compute leaf energy balance
    [LeafState,LeafMassFlux,LeafEnergyFlux] = callEnergyBalance(Constants,CanopyLayer,EnergyOptions,Weather,LeafMassFlux,LeafState);
    LeafState.tLeaf = Error.relax*LeafState.tLeaf + (1-Error.relax)*PreviousLeafState.tLeaf; % Apply relaxation to prevent oscillation
    
    % Compute error in percentage
    Error.gs = abs(PreviousLeafState.gs - LeafState.gs);
    Error.ci = abs(PreviousLeafState.ci - LeafState.ci);
    Error.aNet = abs(PreviousLeafMassFlux.aNet - LeafMassFlux.aNet);
    Error.tLeaf = abs(PreviousLeafState.tLeaf - LeafState.tLeaf);
    
    % Update loop variables
    Error.loop = Error.loop + 1;
    PreviousLeafState = LeafState;
    PreviousLeafMassFlux = LeafMassFlux;
    
%     figure(1)     
%     subplot(221)
%     hold on
%     scatter(Error.loop,LeafMassFlux.aNet)     
%     subplot(222)
%     hold on
%     scatter(Error.loop,LeafState.ci)   
%     subplot(223)
%     hold on
%     scatter(Error.loop,LeafState.gs) 
%     subplot(224)
%     hold on
%     scatter(Error.loop,LeafState.tLeaf)
    
end
if Error.loop >= Error.loopMax
%     disp(strcat("Loop",num2str(loop),"leafLoop",num2str(Error.loop)))
    LeafState.flag = 1;    
else
    LeafState.flag = 0;
end

LeafState.Terror = Error.tLeaf;
LeafState.Aerror = Error.aNet;
LeafState.Cierror = Error.ci;
LeafState.Gserror = Error.gs;
end
