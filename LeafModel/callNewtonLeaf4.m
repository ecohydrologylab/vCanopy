function   [LeafMassFlux,LeafEnergyFlux,LeafState] = callNewtonLeaf4(Constants,...
    LeafBoundaryLayer,EnergyOptions,Photosynthesis,Stomata,Weather,CanopyLayer,...
    LeafState,LeafMassFlux,LeafEnergyFlux,loop)

% Temperature Iteration parameters
Tempiteration = 1;
TempMaxIteration = 100;
CO2MaxIteration = 200;
% MaxNewtonIteration = 200;
lambdaNewton = 0.8;
lambdaTemp = 0.5;
Terror = 10;
Tmaxerror = 0.05;

%% Temperature iteration: Start
while Terror > Tmaxerror && (Tempiteration <= TempMaxIteration)
    
    ComputeTemperatureDependentVariables % Computes temperature depdendent photosynthesis variables
    
    ComputeBoundaryLayerConductance % Computes only the boundary layer conductance
    
    CO2iteration = 1;
    Aerror = 10;
    Cierror = 10;
    Gserror = 10;
    Maxerror = 0.01;
    
    %% CO2 iteration: Start
    while (Aerror > Maxerror || Cierror > Maxerror) && (CO2iteration <= CO2MaxIteration)
        
        %% Does all the preliminary calculation before optimization
        ComputeCO2TransportLimitations % This contains the photosynthesis module and the portion of boundary layer conductance that CO2 depdendent
        
        %% Newton Raphson iteration
%         if CO2iteration <= MaxNewtonIteration
            
            % Scaling for variables
            SX = [10;1000;100;0.1;100];
            % Scaling for functions
            SFx = [10;1000;100;0.1;100];
            
            % Variable Vector
            X = [LeafMassFlux.aNet;LeafState.cbs;LeafState.ci;LeafState.gs;LeafState.cb]./SX;
            
            f1 = @(X) (X(1)*SX(1)) - (1 - LeafState.GammaStar/(X(2)*SX(2)))*...
                ((X(2)*SX(2))*a1/((X(2)*SX(2)) + b1)) + LeafMassFlux.rd;
            f2 = @(X) (X(2)*SX(2)) - (X(3)*SX(3)) - (LeafMassFlux.vp - (X(1)*SX(1))...
                - LeafMassFlux.rm)/Photosynthesis.gbs;
            f3 = @(X) (X(3)*SX(3)) - (X(5)*SX(5)) + 1.6*(X(1)*SX(1))/(X(4)*SX(4));
            if Stomata.model == "BB"
                % Ball-Berry Model
                f4 = @(X) (X(4)*SX(4)) - Stomata.intercept - Stomata.slope*(X(1)*SX(1))*...
                    LeafState.eb/LeafState.ei/((X(5)*SX(5))-LeafState.Gamma);
            elseif Stomata.model == "BM"
                % Belinda Model
                f4 = @(X) (X(4)*SX(4)) - ( Stomata.intercept + (1 + Stomata.slope/(0.001*(LeafState.ei -...
                    LeafState.eb))^0.5)*(X(1)*SX(1))/(X(5)*SX(5)) );
            elseif Stomata.model == "BBL"
                % Ball-Berry-Leunning Model
                f4 = @(X) (X(4)*SX(4)) - Stomata.intercept - Stomata.slope*(X(1)*SX(1))/...
                    (1+(LeafState.ei-LeafState.eb)/(1000*Stomata.D0))/((X(5)*SX(5))-LeafState.Gamma);
            end
            f5 = @(X) (X(5)*SX(5)) - Weather.ca + 1.37*(X(1)*SX(1))/(LeafState.gb*Stomata.Mb);
            
            % Function Vector
            F = @(X) [f1(X);f2(X);f3(X);f4(X);f5(X)]./SFx;
            
            % Stomatal conductance checks
            if LeafState.gs <= Stomata.intercept
                LeafState.gs = Stomata.intercept;
                df4_dX1 = 0;
                df4_dX5 = 0;
            else
                if Stomata.model == "BB"
                    df4_dX1 = -Stomata.slope*(LeafState.eb/LeafState.ei)*(SX(1)/(X(5)*SX(5)-LeafState.Gamma));
                    df4_dX5 = Stomata.slope*(LeafState.eb/LeafState.ei)*(SX(1)*X(1)*SX(5)/(X(5)*SX(5)-LeafState.Gamma)^2);
                elseif Stomata.model == "BM"
                    df4_dX1 = -SX(1)/(X(5)*SX(5))*(1 + Stomata.slope/(0.001*(LeafState.ei -...
                    LeafState.eb))^0.5);
                    df4_dX5 = SX(1)*X(1)/(X(5)^2*SX(5))*(1 + Stomata.slope/(0.001*(LeafState.ei -...
                    LeafState.eb))^0.5);
                elseif Stomata.model == "BBL"
                    df4_dX1 = -Stomata.slope/(1+(LeafState.ei-LeafState.eb)/(1000*Stomata.D0))...
                        *(SX(1)/(X(5)*SX(5)-LeafState.Gamma));
                    df4_dX5 = Stomata.slope/(1+(LeafState.ei-LeafState.eb)/(1000*Stomata.D0))...
                        *(SX(1)*X(1)*SX(5)/(X(5)*SX(5)-LeafState.Gamma)^2);                    
                end
            end
            
            % Derivative of function related to Anet
            df1_dX1 = SX(1);
            df1_dX2 = -(SX(2)*a1*(LeafState.GammaStar + b1)/(b1 + X(2)*SX(2))^2);
            df1_dX3 = 0;
            df1_dX4 = 0;
            df1_dX5 = 0;
            
            % Derivative of function related to Cbs
            df2_dX1 = SX(1)/Photosynthesis.gbs;
            df2_dX2 = SX(2);
            df2_dX3 = -SX(3);
            df2_dX4 = 0;
            df2_dX5 = 0;
            
            % Derivative of function related to Ci
            df3_dX1 = 1.6*SX(1)/(SX(4)*X(4));
            df3_dX2 = 0;
            df3_dX3 = SX(3);
            df3_dX4 = -1.6*SX(1)*X(1)/(SX(4)*X(4)^2);
            df3_dX5 = -SX(5);
            
            % Derivative of function related to gs
            df4_dX2 = 0;
            df4_dX3 = 0;
            df4_dX4 = SX(4);
            
            % Derivative of function related to Cb
            df5_dX1 = 1.37*SX(1)/(LeafState.gb*Stomata.Mb);
            df5_dX2 = 0;
            df5_dX3 = 0;
            df5_dX4 = 0;
            df5_dX5 = SX(5);
            
            % Jacobian Matrix
            dF = [[df1_dX1,df1_dX2,df1_dX3,df1_dX4,df1_dX5]/SFx(1);
                [df2_dX1,df2_dX2,df2_dX3,df2_dX4,df2_dX5]/SFx(2);
                [df3_dX1,df3_dX2,df3_dX3,df3_dX4,df3_dX5]/SFx(3);
                [df4_dX1,df4_dX2,df4_dX3,df4_dX4,df4_dX5]/SFx(4);
                [df5_dX1,df5_dX2,df5_dX3,df5_dX4,df5_dX5]/SFx(5) ];
            
            % Proceeding with the iterations
            try
                if rank(dF) == size(dF,1)
                    X1 = X - dF\F(X);
                else
                    disp(strcat("Singular matrix loop: ",num2str(loop)))
                    X1 = X - pinv(dF)*F(X);
                end
                if norm(F(X1)) > norm(F(X)) %sum(F(X1).^2).^0.5 > sum(F(X).^2).^0.5
                    relax = 0.5;
                    X1 = X - relax*dF\F(X);
                end
            catch
                disp(strcat("Somethings wrong in Input loop: ",num2str(loop)))
            end
            
            if CO2iteration > round(CO2MaxIteration/2)
                lambdaNewton = max(0.1,lambdaNewton - 0.01);
            end
            % Successive under relaxation
            X1 = X1*lambdaNewton + X*(1-lambdaNewton); % Applying relaxation
            
%         elseif CO2iteration > MaxNewtonIteration
%             lambdaNewton = 0.1;
%             % Variable Vector
%             X = [LeafMassFlux.aNet;LeafState.cbs;LeafState.ci;LeafState.gs;LeafState.cb]./SX;
%             
%             %             %% Broyden's method
%             %             if CO2iteration == MaxNewtonIteration + 1
%             %                 B = dF + rand(5).^2; % Initialize the Jacobi matrix
%             %             end
%             %             deltaX = -B\F(X);
%             %             X1 = X + deltaX;
%             %             deltaF = F(X1)-F(X);
%             %             B = B + (deltaF - B*deltaX)*deltaX'/(deltaX'*deltaX);
%             
%             %% Secant method
%             if CO2iteration == MaxNewtonIteration + 1
%                 X1 = X + X/10; % Initialize X1 again with random vars
%             end
%             F1 = F(X1);
%             F0 = F(X);
%             deltaF = (F1 - F0);
%             deltaX = (X1 - X) + 1E-3*ones(5,1);
%             J = deltaF ./ deltaX;
%             X1 = X - lambdaNewton*J;
%             if norm(F(X1)) > norm(F(X))
%                 X1 = X - 0.1*J;
%             end
%             
%             % Successive under relaxation
%             X1 = X1*lambdaNewton + X*(1-lambdaNewton); % Applying relaxation
%         end
        
        if CO2iteration >= CO2MaxIteration
%             F(X1)
        end
        
        % Updating the variables
        Anet = X1(1)*SX(1);
        Cbs = X1(2)*SX(2);
        Ci = X1(3)*SX(3);
        gs = X1(4)*SX(4);
        Cb = X1(5)*SX(5);
        
        % Error calculation
        Aerror = abs(Anet-LeafMassFlux.aNet);
        Cierror = abs(Ci-LeafState.ci);
        Gserror = abs(gs-LeafState.gs);
        
        % Updating the variable
        LeafMassFlux.aNet = Anet; LeafState.cbs = Cbs;
        LeafState.ci = Ci; LeafState.gs = gs; LeafState.cb = Cb;
        
        % Checking the variables
        %% Ball Berry equation for stomatal conductance [mol m-2 sec-1]
        if LeafState.gs < Stomata.intercept || LeafState.gs > 10
            LeafState.gs = Stomata.intercept; % Stomatal conductance for vapour [mol m-2 s-1]
        end
        if LeafState.ci < 0 || LeafState.ci > 2*Weather.ca
            LeafState.ci = 0.7*Weather.ca;
        end
        if LeafState.cb < 0 || LeafState.cb > 2*Weather.ca
            LeafState.cb = 0.8*Weather.ca;
        end
        if LeafState.cbs < 0 || LeafState.cbs > 100*Weather.ca
            LeafState.cbs = 10*LeafState.ci;
        end
        if LeafMassFlux.aNet < -5 || LeafMassFlux.aNet > LeafMassFlux.vc
            LeafMassFlux.aNet = LeafState.gs*(LeafState.cb - LeafState.ci)/1.6;
        end
        
        
        % Storing the error
        LeafState.Aerror = Aerror;
        LeafState.Cierror = Cierror;
        LeafState.Gserror = Gserror;
        
        CO2iteration = CO2iteration + 1;    

    end
    %% CO2 Iteration: End
    
    % Compute leaf temperature
    oldTemp = LeafState.tLeaf;
    [LeafState,LeafMassFlux,LeafEnergyFlux] = callEnergyBalance(Constants,CanopyLayer,EnergyOptions,Weather,LeafMassFlux,LeafState);
    if EnergyOptions.switch == 1
        LeafState.tLeaf = LeafState.tLeaf*lambdaTemp + oldTemp*(1-lambdaTemp);
    end
    Terror = abs(oldTemp-LeafState.tLeaf);
    Tempiteration = Tempiteration + 1;
    
    LeafState.Terror = Terror;
end
%% Temperature iteration: End

if Tempiteration > TempMaxIteration
    LeafState.tLeaf = 0.5*(LeafState.tLeaf + oldTemp);
    LeafState.flag = 1;
%     disp("Leaf Code: Leaf temperature approximated")
end

if CO2iteration > CO2MaxIteration-1
    LeafState.flag = 1;
%     disp(strcat("Leaf Code: Carbon iteration not converged-",num2str(loop)))
else
    LeafState.flag = 0;
end

%% Post processing
% Calculating Leakage
LeafMassFlux.L_bs = Photosynthesis.gbs*(LeafState.cbs-LeafState.ci);
% Photosynthesis (CO2 Fluxes)
LeafMassFlux.aGross = LeafMassFlux.aNet + LeafMassFlux.rd; % Gross rate of CO2 uptake per unit area [u moles m-2 s-1]

end