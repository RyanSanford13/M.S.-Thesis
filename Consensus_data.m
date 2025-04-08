% Read in Consensus Data and other constants needed for Dose Calcs and convert to .mat files for easy reading 

% Length of source
L = 0.35; % cm

% Dose rate constant 
Lambda = 1.117; % cGy / hr-U

% Nominal source strength 
A = 10000; % mCi

% Exposure Rate Constant 
Gamma = 4.69; % R-cm^2 / mCi-hr

% Reference Dose Rate
D_0 = readtable('Consensus Data/192ir-hdr_gammamed_plus.xls', 'Sheet', 1, 'Range', 'AD21:AD21');
D_0 = table2array(D_0);
D_0 = D_0(1);

% Radial Dose Function
gLr = readtable('Consensus Data/192ir-hdr_gammamed_plus.xls', 'Sheet', 1, 'Range', 'B12:C25');
gLr.Properties.VariableNames = {'r (cm)', 'gL(r), L = .35'};

% Anisotropy Function
Frtheta = readtable('Consensus Data/192ir-hdr_gammamed_plus.xls', 'Sheet', 1, 'Range', 'E11:W50');
Frtheta = Frtheta(2:40, 2:19);
Frtheta.Properties.RowNames = {'t = 0', 't = 1', 't = 2', 't = 3', 't = 4', 't = 5', 't = 6', 't = 7', 't = 8', 't = 9', 't = 10', 't = 15', 't = 20', 't = 30', 't = 40', 't = 50', 't = 60', 't = 70', 't = 80', 't = 90', 't = 100', 't = 110', 't = 120', 't = 130', 't = 140', 't = 150', 't = 160', 't = 165', 't = 170', 't = 171', 't = 172', 't = 173', 't = 174', 't = 175', 't = 176', 't = 177', 't = 178', 't = 179', 't = 180'};
Frtheta.Properties.VariableNames = {'r = 0', 'r = 0.2', 'r = 0.4', 'r = 0.6', 'r = 0.8', 'r = 1.0', 'r = 1.25', 'r = 1.5', 'r = 1.75', 'r = 2.0', 'r = 2.5', ' r = 3.0', 'r = 3.5', 'r = 4.0', 'r = 5.0', 'r = 6.0', 'r = 8.0', 'r = 10.0'};

save('L.mat', 'L');
save('Lambda.mat', 'Lambda');
save('A.mat', 'A');
save('Gamma.mat', 'Gamma');
save('D_0.mat', 'D_0');
save('gLr.mat', 'gLr');
save('Frtheta.mat', 'Frtheta');

