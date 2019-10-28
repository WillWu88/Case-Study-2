%SIR Model
xt = [2000; 0; 0; 0];
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
xtpop = [10000; 0; 0; 0; 5000; 0; 0; 0; 15000; 0; 0; 0];

updateMatrixpop = [0.65 0.04 0 0 0    0    0 0 0    0   0 0 ;...
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
                   0    0    0 0 0    0    0 0 0    0.2 0 1];
max_iter = 100;
time_series = 1:max_iter;

susceptible1 = zeros(max_iter, 1);
infected1 = zeros(max_iter, 1);
recovered1 = zeros(max_iter, 1);
dead1 = zeros(max_iter, 1);

xt_temp = num2cell(xtpop(1:4));
[susceptible1(1), infected1(1), recovered1(1), dead1(1)] = deal(xt_temp{:});

for i = 2:max_iter
    xtpop(1:4) = update_xtpop(updateMatrixpop(1:4,1:4), xtpop(1:4));
    xt_temp = num2cell(xtpop(1:4));
    [susceptible1(i), infected1(i), recovered1(i), dead1(i)] = deal(xt_temp{:});
end




susceptible2 = zeros(max_iter, 1);
infected2 = zeros(max_iter, 1);
recovered2 = zeros(max_iter, 1);
dead2 = zeros(max_iter, 1);

xt_temp = num2cell(xtpop(5:8));
[susceptible2(1), infected2(1), recovered2(1), dead2(1)] = deal(xt_temp{:});

for i = 2:max_iter
    xtpop(5:8) = update_xtpop(updateMatrixpop(5:8,5:8), xtpop(5:8));
    xt_temp = num2cell(xtpop(5:8));
    [susceptible2(i), infected2(i), recovered2(i), dead2(i)] = deal(xt_temp{:});
end



susceptible3 = zeros(max_iter, 1);
infected3 = zeros(max_iter, 1);
recovered3 = zeros(max_iter, 1);
dead3 = zeros(max_iter, 1);

xt_temp = num2cell(xtpop(9:12));
[susceptible3(1), infected3(1), recovered3(1), dead3(1)] = deal(xt_temp{:});

for i = 2:max_iter
    xtpop(9:12) = update_xtpop(updateMatrixpop(9:12,9:12), xtpop(9:12));
    xt_temp = num2cell(xtpop(9:12));
    [susceptible3(i), infected3(i), recovered3(i), dead3(i)] = deal(xt_temp{:});
end

hold on;
figure();
subplot(3,1,1);
hold on;
plot(time_series, susceptible1);
plot(time_series, infected1);
plot(time_series, recovered1);
plot(time_series, dead1);
title('SIRD Plot 1 over 500 iterations');
legend('Susceptible','Infected','Recovered','Dead');

subplot(3,1,2);
hold on;
plot(time_series, susceptible2);
plot(time_series, infected2);
plot(time_series, recovered2);
plot(time_series, dead2);
title('SIRD Plot 2 over 500 iterations');
legend('Susceptible','Infected','Recovered','Dead');

subplot(3,1,3);
hold on;
plot(time_series, susceptible3);
plot(time_series, infected3);
plot(time_series, recovered3);
plot(time_series, dead3);
title('SIRD Plot 3 over 500 iterations');
legend('Susceptible','Infected','Recovered','Dead');
hold off;
%%
%Three population SIRD Model
%Traveling populations










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


