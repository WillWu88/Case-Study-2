%SIR Model
xt = [1; 0; 0; 0];
updateMatrix = [0.97 0.04 0 0; 0.03 0.85 0 0; 0 0.1 1 0; 0 0.01 0 1];
max_iter = 200;
time_series = 1:max_iter;

% initializing population groups by city: STL, NYC and Beijing
% as a n by 4 matrix, each column denotes a category such as infected

% setting population as 1 for now, fill in for future reference
STL = zeros(max_iter, 4);
STLPop = 1;
xtSTL = xt * STLPop;

NYC = zeros(max_iter, 4);
NYCPop = 1;
xtNYC = xt * NYCPop;

BJC = zeros(max_iter, 4);
BJCPop = 1;
xtBJC = xt * BJCPop;

umSTL = updateMatrix;
umNYC = [0.85 0.05 0 0; 0.15 0.8 0 0; 0 0.05 1 0; 0 0.1 0 1];
umBJC = [0.6 0.1 0 0; 0.4 0.6 0 0; 0 0.1 1 0; 0 0.2 0 1];

% below
sTravel_rate = [0 0.02 0.01; 0.03 0 0.05; 0.06 0.04 0];
iTravel_rate = [0 0.01 0; 0.08 0 0; 0 0 0];
rTravel_rate = zeros(3, 3);
Travel_rate = cat(3, sTravel_rate, iTravel_rate, rTravel_rate);
% Travel_rate is a 3*3*3 matrix, first layer for s, second for i, third for r
% each entry, Tij, denotes the travelling rate from city i to city j

genUMISO = genMultiUpdateISO(umSTL, umNYC, umBJC);
genUMCON = genMultiUpdateCON(umSTL, umNYC, umBJC, Travel_rate);

% run section ------------------------------------------------------------------
for scenario = 1:4
    xt_constant = [1, 0, 0, 0];

    if scenario == 1
        % first scenario: single population, no vaccine
        xt = xtSTL;
        for i = 1:max_iter
            % simulate single city: STL
            
            xt_temp = num2cell(xt);
            [STL(i, 1), STL(i, 2), STL(i, 3), STL(i, 4)] = deal(xt_temp{:});

            xt = update_xt(updateMatrix, xt);
        end

        figure('Name', 'SIRD Plot: Single Population, No Vaccine');
        hold on;
        plot(time_series, STL(:, 1))
        plot(time_series, STL(:, 2))
        plot(time_series, STL(:, 3))
        plot(time_series, STL(:, 4))
        legend('Susceptible', 'Infected', 'Recovered', 'Dead');
        hold off;

    elseif scenario == 2

        STL = zeros(max_iter, 4);

        for i = 1:max_iter
            %simulate multiple city: STL, NYC and BJC
            if (i == 1)
                % concatonating the 3 state vectors into a 12*1 state vector
                xt = cat(1, xtSTL, xtNYC, xtBJC);
            end

            xt_temp = num2cell(xt);
            [STL(i, 1), STL(i, 2), STL(i, 3), STL(i, 4)] = deal(xt_temp{1:4});
            [NYC(i, 1), NYC(i, 2), NYC(i, 3), NYC(i, 4)] = deal(xt_temp{5:8});
            [BJC(i, 1), BJC(i, 2), BJC(i, 3), BJC(i, 4)] = deal(xt_temp{9:12});

            %isolated scenario
            xt = update_xt(genUMISO, xt);
        end

        plotMultipleSIRD(STL, NYC, BJC, time_series, 'SIRD Plot of STL, NYC and BJC: Isolated');

    elseif scenario == 3

        STL = zeros(max_iter, 4);
        NYC = zeros(max_iter, 4);
        BJC = zeros(max_iter, 4);

        for i = 1:max_iter
            %simulate multiple city: STL, NYC and BJC
            if (i == 1)
                % concatonating the 3 state vectors into a 12*1 state vector
                xt = cat(1, xtSTL, xtNYC, xtBJC);
            end

            xt_temp = num2cell(xt);
            [STL(i, 1), STL(i, 2), STL(i, 3), STL(i, 4)] = deal(xt_temp{1:4});
            [NYC(i, 1), NYC(i, 2), NYC(i, 3), NYC(i, 4)] = deal(xt_temp{5:8});
            [BJC(i, 1), BJC(i, 2), BJC(i, 3), BJC(i, 4)] = deal(xt_temp{9:12});
            % connected scenario
            xt = update_xt(genUMCON, xt);
        end

        plotMultipleSIRD(STL, NYC, BJC, time_series, 'SIRD Plot of STL, NYC and BJC: Connected');

    else
        STL = zeros(max_iter, 4);
        xt = xtSTL;

        for i = 1:max_iter

            xt_temp = num2cell(xt);
            [STL(i, 1), STL(i, 2), STL(i, 3), STL(i, 4)] = deal(xt_temp{:});

            if (i == 20)
                updateMatrix(3, 1) = 0.05;
                updateMatrix(1, 1) = updateMatrix(1, 1) - 0.05;
            end

            xt = update_xt(updateMatrix, xt);
        end

        figure('Name', 'SIRD Plot: Single Population, Vaccinated')
        hold on;
        plot(time_series, STL(:, 1))
        plot(time_series, STL(:, 2))
        plot(time_series, STL(:, 3))
        plot(time_series, STL(:, 4))
        legend('Susceptible', 'Infected', 'Recovered', 'Dead');
        hold off;
    end

end

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
    for layerIndex = 1:l% cycling through three category
        for rindex = 1:r
            for cindex = 1:c
                if rindex == cindex
                    index = 4 * cindex + layerIndex - 4;
                    entry = isoUM(index, index) - sum(trIncidence(rindex, :, layerIndex));
                    isoUM(index, index) = entry;
                else
                    row = 4 * rindex + layerIndex - 4;
                    column = 4 * cindex + layerIndex - 4;
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

function updateMatrix = vaccineIntroMulti(umCity1, umCity2, umCity3, trIncidence, state)
    umCity1(3, 1) = 0.05;
    umCity2(3, 1) = 0.05;
    umCity3(3, 1) = 0.05;

    umCity1(1, 1) = umCity1(1, 1) - 0.05;
    umCity2(1, 1) = umCity2(1, 1) - 0.05;
    umCity3(1, 1) = umCity3(1, 1) - 0.05;

    if (state == 1)
        updateMatrix = genMultiUpdateISO(umCity1, umCity2, umCity3);
    else
        updateMatrix = genMultiUpdateCON(umCity1, umCity2, umCity3, trIncidence);
    end

end

function plotMultipleSIRD(city1, city2, city3, time_series, plot_label)
    figure('Name', plot_label);
    subplot(3, 1, 1);
    hold on;
    plot(time_series, city1(:, 1));
    plot(time_series, city1(:, 2));
    plot(time_series, city1(:, 3));
    plot(time_series, city1(:, 4));
    title('city1');
    legend('Susceptible', 'Infected', 'Recovered', 'Dead');
    hold off;

    subplot(3, 1, 2);
    hold on;
    plot(time_series, city2(:, 1));
    plot(time_series, city2(:, 2));
    plot(time_series, city2(:, 3));
    plot(time_series, city2(:, 4));
    title('city2');
    legend('Susceptible', 'Infected', 'Recovered', 'Dead');
    hold off;

    subplot(3, 1, 3);
    hold on;
    plot(time_series, city3(:, 1));
    plot(time_series, city3(:, 2));
    plot(time_series, city3(:, 3));
    plot(time_series, city3(:, 4));
    title('city3');
    legend('Susceptible', 'Infected', 'Recovered', 'Dead');
    hold off;
end
