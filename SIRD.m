%SIR Model
xt = [1; 0; 0; 0];
updateMatrix = [0.95 0.04 0 0; 0.05 0.85 0 0; 0 0.1 1 0; 0 0.01 0 1];
max_iter = 500;
time_series = 1:max_iter;

% initializing population groups by city: STL, NYC and Beijing
sSTL = zeros(max_iter, 1);
iSTL = zeros(max_iter, 1);
rSTL = zeros(max_iter, 1);
dSTL = zeros(max_iter, 1);
sNYC = zeros(max_iter, 1);
iNYC = zeros(max_iter, 1);
rNYC = zeros(max_iter, 1);
dNYC = zeros(max_iter, 1);
sBJC = zeros(max_iter, 1);
iBJC = zeros(max_iter, 1);
rBJC = zeros(max_iter, 1);
dBJC = zeros(max_iter, 1);

umSTL = updateMatrix;
umNYC = [0.85 0.05 0 0; 0.15 0.8 0 0; 0 0.05 1 0; 0 0.1 0 1];
umBJC = [0.6 0.1 0 0; 0.4 0.6 0 0; 0 0.1 1 0; 0 0.2 0 1];

% below 
travel_rate = [3sy];

genUMISO = genMultiUpdateISO(umSTL, umNYC, umBJC);

xt_temp = num2cell(xt);
[sSTL(1), iSTL(1), rSTL(1), dSTL(1)] = deal(xt_temp{:});
% run section ------------------------------------------------------------------
for i = 2:max_iter
    % simulate single city: STL
    xt = update_xt(updateMatrix, xt);
    xt_temp = num2cell(xt);
    [sSTL(i), iSTL(i), rSTL(i), dSTL(i)] = deal(xt_temp{:});
end

figure('Name', 'SIRD Plot over 500 iterations')
hold on;
plot(time_series, sSTL)
plot(time_series, iSTL)
plot(time_series, rSTL)
plot(time_series, dSTL)
legend('Susceptible', 'Infected', 'Recovered', 'Dead');

function xtPlusOne = update_xt(updateMatrix, xt)
    % takes two arguments: an update matrix and a target matrix
    % returns the next state matrix in time
    xtPlusOne = updateMatrix * xt;
end

function updateMatrix = genMultiUpdateISO(umCity1, umCity2, umCity3)
    % takes the update matrices for the three cities provided
    % given the isolated scenario, generated a 12*12 updated matrix
    % assume all matrices provided are 4*4
    filler = zeros(size(umCity1));
    updateMatrix = cat(1, cat(2, umCity1, filler, filler), cat(2, filler, umCity2, filler), cat(2, filler, filler, umCity3));
end

function updateMatrix = genMultiUpdateCON(umCity1, umCity2, umCity3, trIncidence)
    
end