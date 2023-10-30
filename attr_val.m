%Omolemo Kegakilwe KGKOMO001 EEE4022S Project Model
%computing utilities for users

classdef attr_val
    methods(Static)

        % Cost Attribute Value
        function cost_vals = cost(nwType)
            nwCostValues = struct( 'G3', 0.3843, 'G4', 0.5816, 'G5', (0.7522));
            
            % Initialize an array to store cost values
            cost_vals = zeros(length(nwType),1);

            % Iterate through network types and calculate cost for each
            for i = 1:numel(nwType)
                % Get the cost value for the specified network type
                networkCost = nwCostValues.(nwType{i});
                
                cost_vals(i) = 1-networkCost;
                cost_vals(i) = max(0, min(1, cost_vals(i)));
            end
        end

        % Data Rate Attribute Value
        function rate_vals = datarate(nwType)
            nwRateRanges = struct('G3', [0.5 5], 'G4', [1 50], 'G5', [50 450]);
            rate_vals = zeros(length(nwType),1);
                    
            for i = 1:numel(nwType) 
                nwRange = nwRateRanges.(nwType{i});
                nwRate = normrnd((nwRange(1) + nwRange(2)) / 2, (nwRange(2) - nwRange(1)) / 4);
                rate_vals(i) = max(0,min(1,real(log(nwRate)/log(450))));
            end
        end
        
        function delay_vals = delay(nwType)
            nwDelayRanges = struct('G3', [100 500], 'G4', [5 100], 'G5', [1 4]);

            delay_vals=zeros(length(nwType),1);
            
            for i =1:length(nwType)
                nwRange = nwDelayRanges.(nwType{i});
                nwDelay = normrnd((nwRange(1) + nwRange(2)) / 2, (nwRange(2) - nwRange(1)) / 4);
                delay_vals(i) = max(0,min(1,1-real(log(nwDelay)/log(500))));
            end
        end

        function covrge_vals = coverage(nwType)
            nwcoverageValues = struct('G3', 1.00, 'G4', 0.99, 'G5', 0.30);
            covrge_vals = zeros(length(nwType),1);
            
            for i = 1:numel(nwType)
                 nwcovrge= nwcoverageValues.(nwType{i});
                covrge_vals(i) = max(0,min(1,nwcovrge));
            end
        end

        function attr_vals= all(nwType)
            attr_vals= [attr_val.cost(nwType), attr_val.datarate(nwType),...
            attr_val.delay(nwType), attr_val.coverage(nwType)];
        end

        function [weighted_attr_vals,user_u]=utility(weights, attr_vals)                
            user_u = zeros(length(attr_vals),1);
            weighted_attr_vals= zeros(length(attr_vals),4);
            
            % Calculate the utility for each row
            for i = 1:length(weights)
               weighted_attr_vals(i,:)=weights(i, :) .* attr_vals(i, :);

               user_u(i) = sum(weighted_attr_vals(i,:));
            end
        end

        function [weighted_incentive_attr_vals,incentive_utilites] = incentive_utility(weights, x_attr_vals,...
                g_attr_vals,max_incentive)
            
            incentive_attr_vals=g_attr_vals;
            incentive_utilites=zeros(length(weights),1);
            weighted_incentive_attr_vals=zeros(length(weights),4);

            % Define the maximum incentive (m) and the weights for the attributes
            m = max_incentive;
            
            for i = 1:length(weights)
                u_weights=weights(i,:);
                u_x_attr_vals=x_attr_vals(i,:);
                u_g_attr_vals=g_attr_vals(i,:);

                c=(1-u_g_attr_vals(1));

                alpha=max(0,(u_weights(1)*u_x_attr_vals(1))-(u_weights(1)*u_g_attr_vals(1)))+...
                    max(0,(u_weights(2)*u_x_attr_vals(2))-(u_weights(2)*u_g_attr_vals(2)))+...
                    max(0,(u_weights(3)*u_x_attr_vals(3))-((u_weights(3)*u_g_attr_vals(3))))+...
                    max(0,(u_weights(4)*u_x_attr_vals(4))-((u_weights(4)*u_g_attr_vals(4))));
                
                incentive_attr_vals(i,1)=1-max(c-alpha,c-m);

                weighted_incentive_attr_vals(i,:)=u_weights.* incentive_attr_vals(i, :);

                incentive_utilites(i) = sum(weighted_incentive_attr_vals(i,:));
            end
        end

        function plotUtility(weighted_cost_vals, weighted_delay_vals, weighted_rate_vals,Weighted_covrge_vals, user_utilities, layout_title)
            figure;
            t = tiledlayout('flow');
            
            nexttile
            plot(weighted_cost_vals)
            ylim([0,1])
            title('Weighted Cost Attribute Values')
          
            nexttile
            plot(weighted_delay_vals)
            ylim([0,1])
            title('Weighted Delay Attribute Values')
          
            nexttile
            plot(weighted_rate_vals)
            ylim([0,1])
            title('Weighted Rate Attribute Values')
            
            nexttile
            plot(Weighted_covrge_vals)
            ylim([0,1])
            title('Weighted Coverage Attribute Values')
            
            nexttile([2 4])
            plot(user_utilities)
            ylim([0,1])
            title('Total User Utility')
            
            % Add layout title
            title(t, layout_title)
        end
        
        function result = stats(data)
            minimum = min(data);
            first_quartile = prctile(data, 25);
            median_val = median(data);
            third_quartile = prctile(data, 75);
            maximum = max(data);
            
            std_deviation = std(data);
            mean_val = mean(data);
            
            % Store the computed values in an array
            result = [minimum, first_quartile, median_val, third_quartile, maximum, std_deviation, mean_val];
        end


    end
end
