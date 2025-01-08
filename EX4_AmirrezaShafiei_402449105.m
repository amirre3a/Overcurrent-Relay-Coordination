function EX4_AmirrezaShafiei_402449105
  
    clc;
    clear;

    reportFile = 'GA_RelayCoordination_Report.txt';
    fid = fopen(reportFile, 'w');

    n_relays = 14;
    nvars    = 2 * n_relays;

    TMS_min = 0.05;
    TMS_max = 2.0;

    Ip_min_perRelay = [450; 600; 400; 600; 450; 400; 450; 400; 600; 400; 600; 450; 400; 450];
    Ip_max_perRelay = ones(14,1) * 800;

    lb = zeros(1,nvars);
    ub = zeros(1,nvars);
    for i = 1:n_relays
        lb(i) = TMS_min;
        ub(i) = TMS_max;
        lb(n_relays + i) = Ip_min_perRelay(i);
        ub(n_relays + i) = Ip_max_perRelay(i);
    end

    % Fault Current Data
    If_main = [3232; 5924; 5924; 3556; 3783; 2401; 6109; 6093; 2484; 3883; 3707; 5899; 2991; 5199];
    If_backup = [3232; 996; 1890; 3556; 3783; 2401; 1874; 1890; 2484; 2344; 3707; 987; 2991; 996];
    CTI = 0.2;

    mainBackupPairs = [
         1  6; 2  1; 2  7; 3  2; 4  3; 5  4; 6  5; 6 14; 
         7  5; 7 13; 8  7; 8  9; 9 10; 10 11; 11 12; 12 13; 
        12 14; 13  8; 14  1; 14  9
    ];

    data = struct('n_relays', n_relays, 'If_main', If_main, 'If_backup', If_backup, ...
                  'mainBackupPairs', mainBackupPairs, 'CTI', CTI);

    % GA Options
    opts = optimoptions('ga', ...
        'PopulationSize', 250, ...
        'MaxGenerations', 5000, ...
        'Display', 'off', ...
        'PlotFcn', @gaplotbestf);

    %  GA for Standard Inverse
    fprintf('\n========== GA with Standard Inverse ==========\n');
    fprintf(fid, '\n========== GA with Standard Inverse ==========\n');
    data.a = 0.14;
    data.b = 0.02;
    [x_opt_SI, fval_SI] = ga(@(x) objectiveRelay(x, data), nvars, [], [], [], [], lb, ub, [], opts);

    TMS_SI = x_opt_SI(1:n_relays);
    Ip_SI  = x_opt_SI(n_relays+1:end);
    avg_TMS_SI = mean(TMS_SI);
    avg_Ip_SI = mean(Ip_SI);

    fprintf('\n** Standard Inverse Results **\n');
    fprintf('Best Objective = %.4f\n', fval_SI);
    fprintf('Optimal TMS:\n'); disp(TMS_SI);
    fprintf('Optimal Ip:\n');  disp(Ip_SI);
    fprintf('Average TMS = %.4f\n', avg_TMS_SI);
    fprintf('Average Ip = %.4f\n', avg_Ip_SI);

    fprintf(fid, '** Standard Inverse Results **\n');
    fprintf(fid, 'Best Objective = %.4f\n', fval_SI);
    fprintf(fid, 'Optimal TMS:\n');
    fprintf(fid, '%f\n', TMS_SI);
    fprintf(fid, 'Optimal Ip:\n');
    fprintf(fid, '%f\n', Ip_SI);
    fprintf(fid, 'Average TMS = %.4f\n', avg_TMS_SI);
    fprintf(fid, 'Average Ip = %.4f\n\n', avg_Ip_SI);

    %  GA for Very Inverse
    fprintf('\n========== GA with Very Inverse ==========\n');
    fprintf(fid, '\n========== GA with Very Inverse ==========\n');
    data.a = 13.5;
    data.b = 1.0;
    [x_opt_VI, fval_VI] = ga(@(x) objectiveRelay(x, data), nvars, [], [], [], [], lb, ub, [], opts);

    TMS_VI = x_opt_VI(1:n_relays);
    Ip_VI  = x_opt_VI(n_relays+1:end);
    avg_TMS_VI = mean(TMS_VI);
    avg_Ip_VI = mean(Ip_VI);

    fprintf('\n** Very Inverse Results **\n');
    fprintf('Best Objective = %.4f\n', fval_VI);
    fprintf('Optimal TMS:\n'); disp(TMS_VI);
    fprintf('Optimal Ip:\n');  disp(Ip_VI);
    fprintf('Average TMS = %.4f\n', avg_TMS_VI);
    fprintf('Average Ip = %.4f\n', avg_Ip_VI);

    fprintf(fid, '** Very Inverse Results **\n');
    fprintf(fid, 'Best Objective = %.4f\n', fval_VI);
    fprintf(fid, 'Optimal TMS:\n');
    fprintf(fid, '%f\n', TMS_VI);
    fprintf(fid, 'Optimal Ip:\n');
    fprintf(fid, '%f\n', Ip_VI);
    fprintf(fid, 'Average TMS = %.4f\n', avg_TMS_VI);
    fprintf(fid, 'Average Ip = %.4f\n\n', avg_Ip_VI);

    fclose(fid);
    fprintf('Results saved to %s\n', reportFile);
end

function OF = objectiveRelay(x, data)
    n_relays = data.n_relays;
    If_main  = data.If_main;
    If_backup = data.If_backup;
    pairs = data.mainBackupPairs;
    a = data.a;
    b = data.b;
    CTI = data.CTI;

    TMS = x(1:n_relays);
    Ip = x(n_relays+1:end);

    t_main = zeros(n_relays, 1);
    for i = 1:n_relays
        ratio = (If_main(i) / Ip(i))^b;
        if ratio <= 1
            t_main(i) = 9999;
        else
            t_main(i) = TMS(i) * (a / (ratio - 1));
        end
    end

    T_total = sum(t_main);

    bigPenalty = 100;
    penaltyCount = 0;
    for k = 1:size(pairs, 1)
        mainR = pairs(k, 1);
        backR = pairs(k, 2);

        ratio_b = (If_backup(backR) / Ip(backR))^b;
        if ratio_b <= 1
            t_back = 9999;
        else
            t_back = TMS(backR) * (a / (ratio_b - 1));
        end

        if t_back - t_main(mainR) < CTI
            penaltyCount = penaltyCount + 1;
        end
    end

    OF = T_total + bigPenalty * penaltyCount;
end
