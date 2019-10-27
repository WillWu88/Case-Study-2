%SIR Model
xt = [1; 0; 0; 0];
updateMatrix = [0.95 0.04 0 0;...
                0.05 0.85 0 0;...
                0    0.1  1 0;...
                0    0.01 0 1];
max_iter = 500;
time_series = 1:max_iter;

susceptible = zeros(max_iter, 1);
infected = zeros(max_iter, 1);
recovered = zeros(max_iter, 1);
dead = zeros(max_iter, 1);

xt_temp = num2cell(xt);
[susceptible(1), infected(1), recovered(1), dead(1)] = deal(xt_temp{:});

for i = 2:max_iter
    xt = update_xt(updateMatrix, xt);
    xt_temp = num2cell(xt);
    [susceptible(i), infected(i), recovered(i), dead(i)] = deal(xt_temp{:});
end


figure('Name', 'SIRD Plot over 500 iterations')
hold on;
plot(time_series, susceptible);
plot(time_series, infected);
plot(time_series, recovered);
plot(time_series, dead);
legend('Susceptible','Infected','Recovered','Dead');



%%
%Three population SIRD Model
%Non-traveling populations
xtpop = [1; 0; 0; 0; 1; 0; 0; 0; 1; 0; 0; 0];
updateMatrixpop = [0.95 0.04 0 0 0    0    0 0 0    0   0 0 ;...
                   0.05 0.85 0 0 0    0    0 0 0    0   0 0 ;...
                   0    0.1  1 0 0    0    0 0 0    0   0 0 ;...
                   0    0.01 0 1 0    0    0 0 0    0   0 0 ;...
                   0    0    0 0 0.85 0.05 0 0 0    0   0 0 ;...
                   0    0    0 0 0.15 0.8  0 0 0    0   0 0 ;...
                   0    0    0 0 0    0.05 1 0 0    0   0 0 ;...
                   0    0    0 0 0    0.1  0 1 0    0   0 0 ;...
                   0    0    0 0 0    0    0 0 0.6  0.1 0 0 ;...
                   0    0    0 0 0    0    0 0 0.4  0.6 0 0 ;...
                   0    0    0 0 0    0    0 0 0    0.1 1 0 ;...
                   0    0    0 0 0    0    0 0 0    0.2 0 1];...
max_iter = 500;
time_series = 1:max_iter;

susceptible = zeros(max_iter, 1);
infected = zeros(max_iter, 1);
recovered = zeros(max_iter, 1);
dead = zeros(max_iter, 1);

xt_temp = num2cell(xt);
[susceptible(1), infected(1), recovered(1), dead(1)] = deal(xt_temp{:});

for i = 2:max_iter
    xtpop = update_xtpop(updateMatrixpop, xtpop);
    xt_temp = num2cell(xtpop);
    [susceptible(i), infected(i), recovered(i), dead(i)] = deal(xt_temp{:});
end


figure('Name', 'SIRD Plot over 500 iterations')
hold on;
plot(time_series, susceptible);
plot(time_series, infected);
plot(time_series, recovered);
plot(time_series, dead);
legend('Susceptible','Infected','Recovered','Dead');



function xtPlusOne = update_xt(updateMatrix, xt)
% takes two arguments: an update matrix and a target matrix
% returns the next state matrix in time
    xtPlusOne = updateMatrix * xt;
end

function xtPlusOnepop = update_xtpop(updateMatrixpop, xtpop)
% takes two arguments: an update matrix and a target matrix
% returns the next state matrix in time
    xtPlusOnepop = updateMatrixpop * xtpop;
end

