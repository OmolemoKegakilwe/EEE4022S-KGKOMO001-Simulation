%Omolemo Kegakilwe KGKOMO001 EEE4022S Project Model
%Find different HWN distributions


filepath="figures/results/";
    if ~isfolder(filepath)
        mkdir(filepath);  
    end

numUsers = 10000;

% Defining User Weights/preferances
prob_DLM1 = 0.1; 
prob_DLM2 = 0.2; 
prob_DLM3 = 0.3; 
prob_DLM4 = 0.2;
prob_DLM5 = 0.1;    

categoryLabels = discretize(rand(1, numUsers), [0, prob_DLM1, prob_DLM1 + ...
    prob_DLM2, prob_DLM1 + prob_DLM2 + prob_DLM3, prob_DLM1 + prob_DLM2 +...
    prob_DLM3 + prob_DLM4, 1], {'DLM1', 'DLM2', 'DLM3', 'DLM4', 'DLM5'});
categoryLabels = sort(categoryLabels);

weightsMatrix = zeros(numUsers, 4); 

for i = 1:numUsers
    switch categoryLabels{i}
        case 'DLM1'
            weightsMatrix(i, :) = [1 * rand(), 0.0 + 0.4 * rand(),...
                                   0.0 + 0.4 * rand(), 0.0 + 0.4 * rand()];
        case 'DLM2'
            weightsMatrix(i, :) = [rand(), 0.2 + 0.3 * rand(),...
                                   0.3 + 0.4 * rand(), 0.3 + 0.4 * rand()];
        case 'DLM3'
            weightsMatrix(i, :) = [rand(), 0.3 + 0.4 * rand(),...
                                   0.4 + 0.4 * rand(), 0.4 + 0.4 * rand()];
        case 'DLM4'
            weightsMatrix(i, :) = [rand(), 0.5 + 0.3 * rand(),...
                                   0.5 + 0.3 * rand(), 0.5 + 0.3 * rand()];
        case 'DLM5'
            weightsMatrix(i, :) = [rand(), 0.7 + 0.3 * rand(),...
                                   0.7 + 0.3 * rand(), 0.7  + 0.3 * rand()];
    end
end
normalizedWeights = weightsMatrix ./ sum(weightsMatrix, 2);

%determine baseline values(control measurement)
scale_users=numUsers/(numUsers-numUsers*0.017);
num3G = round((0.58 * numUsers)*scale_users);
num4G = round((0.4 * numUsers)*scale_users);
num5G = round((0.003 * numUsers)*scale_users);

nw_users = [num3G, num4G, num5G];
baseline = [repmat("G3", round(num3G), 1);
                  repmat("G4", num4G, 1);
                  repmat("G5", num5G, 1)];

%calculate and plot baseline total user utility
% Determine attribute values
cost_vals= attr_val.cost(baseline);
rate_vals= attr_val.datarate(baseline);
delay_vals= attr_val.delay(baseline);
covrge_vals=attr_val.coverage(baseline);

baseline_attr_vals = [cost_vals, delay_vals, rate_vals, covrge_vals];
[weighted_baseline_attr_vals,baseline_utilities]=...
    attr_val.utility(normalizedWeights,baseline_attr_vals);
baseline_utilities_stats=attr_val.stats(baseline_utilities);
average_utility = cellfun(@mean, mat2cell(baseline_utilities, nw_users, 1));
counts_baseline = countcats(categorical(baseline, ["G3", "G4", "G5"]));


% user utility on 5G current sub
nw_sub=repmat("G5", numUsers, 1);
attr_vals_5G= attr_val.all(nw_sub);
[weighted_5G_attr_vals,utilities_on_5G]=...
    attr_val.utility(normalizedWeights,attr_vals_5G);
utilities_on_5G_stats=attr_val.stats(baseline_utilities);


utilities={baseline_utilities,utilities_on_5G};
utilities_labels = {'Baseline', '5G'};

% Define an array of reserve prices
max_incentives = [0,0.05,0.1, 0.15, 0.2, 0.25,0.30];

% HWN user distribution with differing threshold values
threshold_levels = [0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, 0.25,0.27,0.30];

% Create a cell array to store results for each reserve price
results = cell(length(max_incentives), 1);


% Iterate through reserve prices
for m_incentive_index = 1:length(max_incentives)
    max_incentive = max_incentives(m_incentive_index);

    % Evaluate utilities with incentive
    [weighted_incentive_attr_vals, utilities_w_5G_incentive] = ...
        attr_val.incentive_utility(normalizedWeights, baseline_attr_vals,...
        attr_vals_5G, max_incentive);
    utilities{end+1}=utilities_w_5G_incentive;
    utilities_labels{end+1}= strcat('5G(max incentive: ',num2str(max_incentive),')');


    % Initialize cell array for this reserve price
    HWN_distributions = cell(size(threshold_levels));

    % Iterate through threshold levels
    for t_level = 1:length(threshold_levels)
        t = threshold_levels(t_level);
        HWN_distribution = baseline;

        % Iterate through each user
        for userIndex = 1:numUsers
            subscription = char(baseline(userIndex));
            migration_nw = str2double(subscription(2));

            if (migration_nw < 5)
                u_i = utilities_w_5G_incentive(userIndex);
                u_b = baseline_utilities(userIndex);
                u_c = utilities_on_5G(userIndex);

                % Check if migration to 5G is justified
                if (u_i > u_b) && (u_i > (t + u_c))
                    migration_nw = 5;
                end
            end

            HWN_distribution(userIndex) = strcat("G", num2str(migration_nw));
        end

        HWN_distributions{t_level} = HWN_distribution;
    end

    % Store the results for this reserve price
    results{m_incentive_index} = HWN_distributions;
end

%plot box-whisker for utilities
filename="Box Plot of User Utilities";
figure;
boxplot(cell2mat(utilities));
title('Box Plot of User Utilities on various Networks');
ylim([0,1]);
xticklabels(utilities_labels);
saveas(gcf, fullfile(filepath, strcat(filename, '.png')));
clear gcf;

allRATcounts=[];
% %distributions for each threshold, reserve price combination
for max_incentiveIdx=1:length(max_incentives)
    HWN_distributions=results{max_incentiveIdx};
    filepath='figures/results/HWNdistributions/';
    
    if ~isfolder(filepath)
        mkdir(filepath);  
    end

    RATcounts=[];
    for HWNIdx=1:length(threshold_levels)
        HWN_distribution=HWN_distributions{HWNIdx};
        RATcounts=[RATcounts;(countcats(categorical(HWN_distribution,...
            ["G3", "G4", "G5"])))'];
        allRATcounts=[allRATcounts;RATcounts];


    end

    filename=strcat('HWN Distribution of Subscriptions with maxincentive(',...
        num2str(max_incentives(max_incentiveIdx)),')');
    
    figure;
    bar(threshold_levels,RATcounts);
    xlabel('Migration Threshold Values');
    ylabel('Number of Cellular Network subscriptions');
    title(strcat('Network Subscriptions as Threshold changes w maxincentive=',...
        num2str(max_incentives(max_incentiveIdx))));
    legend({'3G Subscriptions', '4G Subscriptions', '5G Subscriptions'});
    h = get(gca, 'Children');
    set(h(1), 'FaceColor', 'c');
    set(h(2), 'FaceColor', 'b');
    set(h(3), 'FaceColor', 'g');

    saveas(gcf, fullfile(filepath, strcat(filename, '.png')));
    clear gcf;
end


%Admissible States
%determine baseline values(control measurement)
num4G = round((0.4 * numUsers)*scale_users);
num5G = round((0.003 * numUsers)*scale_users);
num3G = round((0.58 * numUsers)*scale_users);


baseline = [repmat("G3", num3G, 1);repmat("G4", num4G, 1);repmat("G5", num5G, 1)];
% Network Performance Analysis
arrival_rates = [0,2,4,6,8,10,12,14,16,18,20];


%Capacity, Threshold parameters
nw_param_scale=2;
C1=30/nw_param_scale;   C2=60/nw_param_scale;  C3=100/nw_param_scale;
T1=ceil(C1*0.65);      T2=ceil(C2*0.65);      T3=ceil(C3*0.65);

%Basic bandwidth units
b1=4;       b2=2;       b3=1;

max_m = [ceil(T1/b1),ceil(C1/b1),...
    ceil(T2/b2),ceil(C2/b2),...
    ceil(T3/b3),ceil(C3/b3)];
% Initialize a variable to store the admissible states
admss_states1=[];
admss_states2=[];
admss_states3=[];

for m11 = 0:max_m(1)
    for m21 = 0:max_m(2)
        if ((m11+m21)*b1 <= C1) && (m11*b1<=T1)
            admss_states1=[admss_states1;[m11,m21]];
        end
    end
end
for m12 = 0:max_m(3)
    for m22 = 0:max_m(4)
        if ((m12+m22)*b2 <= C2) && (m12*b2<=T2)
            admss_states2=[admss_states2;[m12,m22]];
        end
    end
end

% Get the sizes of the input matrices
size1 = size(admss_states1, 1);
size2 = size(admss_states2, 1);

% Create a matrix with the appropriate dimensions
admss_states12 = zeros(size1 * size2, 4);

% Populate the columns of admss_states12
for i = 1:size1
    admss_states12((i - 1) * size2 + 1 : i * size2, 1) = admss_states1(i, 1);
    admss_states12((i - 1) * size2 + 1 : i * size2, 2) = admss_states1(i, 2);
end

for j = 1:size2
    admss_states12(j:size2:end, 3) = admss_states2(j, 1);
    admss_states12(j:size2:end, 4) = admss_states2(j, 2);
end

clear admss_states1 admss_states2

for m13 = 0:max_m(5)
    for m23 = 0:max_m(6)
        if ((m13+m23)*b3 <= C3) && (m13*b3<=T3)
            admss_states3=[admss_states3;[m13,m23]];
        end
    end
end

% Get the sizes of the input matrices
size12 = size(admss_states12, 1);
size3 = size(admss_states3, 1);

% Create a matrix with the appropriate dimensions
admss_states = zeros(size12 * size3, 6);

% Populate the columns of admss_states123
for i = 1:size12
    admss_states((i - 1) * size3 + 1 : i * size3,1) = admss_states12(i, 1);
    admss_states((i - 1) * size3 + 1 : i * size3,2) = admss_states12(i, 2);
    admss_states((i - 1) * size3 + 1 : i * size3,3) = admss_states12(i, 3);
    admss_states((i - 1) * size3 + 1 : i * size3,4) = admss_states12(i, 4);
end

for j = 1:size3
    admss_states(j:size3:end, 5) = admss_states3(j, 1);
    admss_states(j:size3:end, 6) = admss_states3(j, 2);
end

% Clear admss_states12 from memory
clear admss_states12 admss_states3

for max_incentiveIdx=1:length(max_incentives)
    HWN_distributions=results{max_incentiveIdx};
    filepath='figures/results/NetworkPerformance/';

    utilization_values={};
    labels={};

    call_dropping={};
    call_blocking={};

        
    if ~isfolder(filepath)
        mkdir(filepath);  
    end

    for HWNIdx=1:length(threshold_levels)
        HWN_distribution=HWN_distributions{HWNIdx};
        RATcount=countcats(categorical(HWN_distribution, ["G3", "G4", "G5"]));

        HWN_utilization_values = zeros(length(arrival_rates),1);
        HWN_call_blocking = zeros(length(arrival_rates),1);
        HWN_call_dropping = zeros(length(arrival_rates),1);  

        for arrIdx = 1:length(arrival_rates)
        
            total_xm = arrival_rates(arrIdx); % Set the current total arrival rate
        
           
            % Calculate the arrival rate for each RAT based on proportions
            lambda_arr = zeros(length(RATcount),1);
        
            lambda_arr(1)=total_xm*((RATcount(1)/numUsers)+...
                ((RATcount(2)/numUsers)*(1/2))+...
                ((RATcount(3)/numUsers)*(1/3)));
        
            lambda_arr(2)=total_xm*(((RATcount(2)/numUsers)*(1/2))+...
                ((RATcount(3)/numUsers)*(1/3)));

            lambda_arr(3)=total_xm*((RATcount(3)/numUsers)*(1/3));
        
         
            %service rates
            mu_dep = zeros(length(RATcount),1);

            mu_dep(1)=0.5;
            mu_dep(2)=0.5;
            mu_dep(3)=0.5;

        
            rat_params = [
            struct('xm', lambda_arr(1), 'xn', lambda_arr(1)/5, 'um', mu_dep(1), 'un', mu_dep(1));
            struct('xm', lambda_arr(2), 'xn', lambda_arr(2)/5, 'um', mu_dep(2), 'un', mu_dep(2));
            struct('xm', lambda_arr(3), 'xn', lambda_arr(3)/5, 'um', mu_dep(3), 'un', mu_dep(3));
        ];
            
            
            %poission processes for determining arrival rates
            rho11=rat_params(1).xm/rat_params(1).um;
            rho21=rat_params(1).xn/rat_params(1).un;
            rho12=rat_params(2).xm/rat_params(2).um;
            rho22=rat_params(2).xn/rat_params(2).un;
            rho13=rat_params(3).xm/rat_params(3).um;
            rho23=rat_params(3).xn/rat_params(3).un;
    
           
           rho=[rho11,rho21,rho12,rho22,rho13,rho23]; 
        
            % Calculate the normalization constant G and utilisation in Admissable States
            G = 0;SB=0;SD=0;cntD=0;cntB=0;
            U_s=zeros(length(admss_states),1);
            P_s = zeros(length(admss_states),1);
            for i = 1:size(admss_states, 1)
                m = admss_states(i, :);
                product = (((rho(1)^m(1))*(rho(2)^m(2)))/(factorial(m(1))*factorial(m(2))))*...
                    (((rho(3)^m(3))*(rho(4)^m(4)))/(factorial(m(3))*factorial(m(4))))*...
                    (((rho(5)^m(5))*(rho(6)^m(6)))/(factorial(m(5))*factorial(m(6))));
        
                G = G + product;
                P_s(i) = product;
                U_s(i)=((m(1)+m(2))*b1)+((m(3)+m(4))*b2)+((m(5)+m(6))*b3);
        
                %blocking probability
                if ((b1+(b1*(m(1)+m(2))))>C1) || (b1+(b1*m(1))>T1) ||...
                        (b2+(b2*(m(3)+m(4)))>C2) || (b2+(b2*m(3))>T2)||...
                        (b3+(b3*(m(5)+m(6)))>C3) || (b3+(b3*m(5))>T3)
                
                    SB=SB+product;
                    cntB=cntB+1;
                end
        
                %Dropping Probability
                if (b1+(b1*(m(1)+m(2)))>C1) ||...
                        (b2+(b2*(m(3)+m(4)))>C2)||...
                        (b3+(b3*(m(5)+m(6)))>C3)
                    SD=SD+product;
                end
            end

            HWN_call_blocking(arrIdx)=SB/G;
            HWN_call_dropping(arrIdx)=SD/G;

            P_s=P_s/G;
            U=sum(P_s.*U_s);
            HWN_utilization_values(arrIdx) =U/(C1+C2+C3);
        end

    utilization_values{end+1}=HWN_utilization_values;
    labels{end+1}=strcat('migration thr:',num2str(threshold_levels(HWNIdx)));
    call_blocking{end+1}=HWN_call_blocking;
    call_dropping{end+1}=HWN_call_dropping;

    end    

    filename=strcat('Call Dropping with maxincentive(',num2str(max_incentives(max_incentiveIdx)),')');
    figure;
    plot(arrival_rates,cell2mat(call_dropping));
    ylim([0,1]);
    xlabel('Arrival Rate');
    ylabel('Call Dropping Probability');
    title(strcat('Dropping Probability as Threshold changes w maxincentive=',num2str(max_incentives(max_incentiveIdx))));
    legend(labels);
    saveas(gcf, fullfile(filepath, strcat(filename, '.png')));
    clear gcf;

    filename=strcat('Call Blocking with maxincentive(',num2str(max_incentives(max_incentiveIdx)),')');
    figure;
    plot(arrival_rates,cell2mat(call_blocking));
    ylim([0,1]);
    xlabel('Arrival Rate');
    ylabel('Call Blocking Probability');
    title(strcat('Blocking Probability as Threshold changes w maxincentive=',num2str(max_incentives(max_incentiveIdx))));
    legend(labels);
    saveas(gcf, fullfile(filepath, strcat(filename, '.png')));
    clear gcf;
    
    filename=strcat('Resource Utilisation with maxincentive(',num2str(max_incentives(max_incentiveIdx)),')');
    figure;
    plot(arrival_rates,cell2mat(utilization_values));
    ylim([0,1]);
    xlabel('Arrival Rate');
    ylabel('Resource Utilisation');
    title(strcat('Resource Utilisation as Threshold changes w maxincentive=',num2str(max_incentives(max_incentiveIdx))));
    legend(labels);
    saveas(gcf, fullfile(filepath, strcat(filename, '.png')));
    clear gcf;

end