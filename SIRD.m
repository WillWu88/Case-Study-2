
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
sTravel_rate = [0 0.02 0.01; 0.03 0 0.05; 0.06 0.04 0];
iTravel_rate = [0 0.01 0; 0.08 0 0; 0 0 0];
rTravel_rate = [0 0.01 0.01; 0.01 0 0.01; 0.01 0.01 0];
Travel_rate = cat(3, sTravel_rate, iTravel_rate, rTravel_rate);
% Travel_rate is a 3*3*3 matrix, first layer for s, second for i, third for r
% each entry, Tij, denotes the travelling rate from city i to city j

genUMISO = genMultiUpdateISO(umSTL, umNYC, umBJC);
genUMCON = genMultiUpdateCON(umSTL, umNYC, umBJC, Travel_rate);

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

hold on;
figure();
subplot(3,1,1);
hold on;
plot(time_series, sSTL);
plot(time_series, iSTL);
plot(time_series, rSTL);
plot(time_series, dSTL);
title('SIRD Plot 1 over 500 iterations');
legend('Susceptible','Infected','Recovered','Dead');

subplot(3,1,2);
hold on;
plot(time_series, sNYC);
plot(time_series, iNYC);
plot(time_series, rNYC);
plot(time_series, dNYC);
title('SIRD Plot 2 over 500 iterations');
legend('Susceptible','Infected','Recovered','Dead');

subplot(3,1,3);
hold on;
plot(time_series, sBJC);
plot(time_series, iBJC);
plot(time_series, rBJC);
plot(time_series, dBJC);
title('SIRD Plot 3 over 500 iterations');
legend('Susceptible','Infected','Recovered','Dead');
hold off;

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
    % takes the update matrices for the three cities and the travel rate matrix
    % given that travel is allowed, generate a 12*12 update matrix
    % assume trIncidence is 3*3*3, first layer for s, second for i, third for r
    isoUM = genMultiUpdateISO(umCity1, umCity2, umCity3);
    [r, c, l] = size(trIncidence);
    for layerIndex = 1:l % cycling through three category
        for rindex = 1:r
           for cindex = 1:c
            if rindex == cindex
                index = 4*cindex+layerIndex-4;
                entry = isoUM(index, index) - sum(trIncidence(rindex,:,layerIndex));
                isoUM(index, index) = entry;
            else
                row = 4*rindex+layerIndex-4;
                column = 4*cindex+layerIndex-4;
                isoUM(row, column) = isoUM(row, column) + trIncidence(cindex, rindex, layerIndex);
                % flipped cindex and rindex, since the loop is looping through 
                % rows then columns, but the filling order should be columns
                % then rows
                
            end
           end
        end
    end
    updateMatrix = isoUM;
end
