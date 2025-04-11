function    [Canopy] = callCanopyGrid(Canopy,Weather,Options)

Canopy.deltaZ = (Canopy.hN-Canopy.h0)/(Options.canopyLayers-1); % Discretization [m]
Canopy.zh = (Canopy.h0:Canopy.deltaZ:Canopy.hN).*ones(height(Weather),1); % Canopy layer bounadaries 
Canopy.zc = Canopy.zh - Canopy.deltaZ/2; % Only the centeral nodes
Canopy.zc(:,1) = [];
Canopy.z = [zeros(size(Canopy.zc,1),1) Canopy.zc]; % Central with soil layer
Canopy.zn = [Canopy.z Canopy.zh(1,end)*ones(size(Canopy.z,1),1)]; % Central with soil and weather
% zeros(size(Canopy.z,1),1);
% Canopy.z = [Canopy.z Canopy.h(1,end)*ones(size(Canopy.z,1),1)];
Canopy.deltaZ = zeros(size(Canopy.z));
Canopy.deltaZ(:,1:end) = diff(Canopy.zn')';
% Data1 = table2array(readtable('./Input/LeafOutput.xlsx'));
% LAIz = Data1(1:canopyLayers,6);
% LAIz(1) = 0;
% Canopy.zn = [Canopy.z

% % % for i=1:1:height(Weather)
% % %     Canopy.LAIProfile(i,:) = LAI;
% % %     end
% % % Canopy.LADProfile = Canopy.LAIProfile./Canopy.deltaZ;
% % LAIz = Data1(1:size(Canopy.z,2),6);
% % 
% %     for i=1:1:length(Canopy.z)
% %     Canopy.LAIProfile(i,:) = LAIz';
% %     end
% % Canopy.LADProfile = Canopy.LAIProfile./Canopy.deltaZ;

end