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
xt_pop = [1; 0; 0; 0; 1; 0; 0; 0; 1; 0; 0; 0];
updateMatrix_pop =[0.95 0.04 0 0 0    0    0 0 0    0   0 0 ;...
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

susceptible1 = zeros(max_iter, 1);
infected1 = zeros(max_iter, 1);
recovered1 = zeros(max_iter, 1);
dead1 = zeros(max_iter, 1);

susceptible2 = zeros(max_iter, 1);
infected2 = zeros(max_iter, 1);
recovered2 = zeros(max_iter, 1);
dead2 = zeros(max_iter, 1);

susceptible3 = zeros(max_iter, 1);
infected3 = zeros(max_iter, 1);
recovered3 = zeros(max_iter, 1);
dead3 = zeros(max_iter, 1);

xt_temp = num2cell(xt_pop);
[susceptible1(1), infected1(1), recovered1(1), dead1(1), susceptible2(1), infected2(1), recovered2(1), dead2(1), susceptible3(1), infected3(1), recovered3(1), dead3(1)] = deal(xt_temp{:});

for i = 1:max_iter
    xt_pop = update_xt_pop(updateMatrix_pop, xt_pop);
    xt_temp = num2cell(xt_pop);
    [susceptible1(i), infected1(i), recovered1(i), dead1(i)] = deal(xt_temp{:});
end
for i = 1:max_iter
    xt_pop = update_xt_pop(updateMatrix_pop, xt_pop);
    xt_temp = num2cell(xt_pop);
    [susceptible2(i), infected2(i), recovered2(i), dead2(i)] = deal(xt_temp{:});
end
for i = 1:max_iter
    xt_pop = update_xt_pop(updateMatrix_pop, xt_pop);
    xt_temp = num2cell(xt_pop);
    [susceptible3(i), infected3(i), recovered3(i), dead3(i)] = deal(xt_temp{:});
end

%figure('Name', 'SIRD Plot over 500 iterations')
%hold on;
%plot(time_series, susceptible);
%plot(time_series, infected);
%plot(time_series, recovered);
%plot(time_series, dead);
%legend('Susceptible','Infected','Recovered','Dead');



function xtPlusOne = update_xt(updateMatrix, xt)
% takes two arguments: an update matrix and a target matrix
% returns the next state matrix in time
    xtPlusOne = updateMatrix * xt;
end

function xtPlusOne_pop = update_xt_pop(updateMatrix_pop, xt_pop)
% takes two arguments: an update matrix and a target matrix
% returns the next state matrix in time
    xtPlusOnepop = updateMatrixpop * xtpop;
end

