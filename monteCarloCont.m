function mc = monteCarloCont( samples, priors, numArms, numTypes, allUsers )
%MONTECARLO Summary of this function goes here
%   Detailed explanation goes here

%% Sample beta distributions and fit single beta distribution 
    %Number of samples:
    sampleN = samples ;

    % Collection of the samples across arm and context
    
    sampled_fit = zeros(sampleN, numArms, ceil(numTypes/8)) ;
    alpha = zeros(2,numArms);
    userCount = [(1:numTypes).',histc(allUsers,1:numTypes)];
    topUsers = sortrows(userCount,-2);
    countJ = 1;
    list = topUsers(1:ceil(numTypes/8),1);

    % Amalgamated samples across arms
    sampled_prob = zeros(sampleN,numArms) ;
    for I = 1 : numArms
       for J = list.'
           sampled_fit(:,I,countJ) = random('beta',priors(1,I,J),priors(2,I,J),sampleN,1) ; 
           countJ = countJ + 1;
       end
       tmp = reshape(sampled_fit(:,I,:),[sampleN*ceil(numTypes/8),1]);
       alpha(:,I) = betafit(tmp,0.75);
       countJ = 1;
    end
    
    mc = monteCarloUni(sampleN,alpha,numArms);
    
%% Numerical addition of beta distributions by occurance weight
%     sampled_fit = zeros(sampleN, numArms, ceil(numTypes/4)) ;
%     alpha = zeros(2,numArms);
%     X = 0:.001:0.999;
%     plotCo = zeros(1,length(X));
%     counts = histc(allUsers,1:numTypes);
%     weights = counts / sum(counts);
%     data = [];
%     
%     % Amalgamated across userTypes based on weighted occurrence
%     for I = 1 : numArms
%        for J = 1 : numTypes
%            plotCo = plotCo + weights(J)*betapdf(X,priors(1,I,J),priors(2,I,J));
%        end
%        plotCo = floor(plotCo);
%        for J = 1 : length(plotCo)
%            if plotCo(J) ~= 0
%                 tmp = 0.001*J*ones(1,plotCo(J));
%                 data = [data tmp];
%            end
%        end
%        alpha(:,I) = ceil(betafit(data,0.1));
%     end
%     
%     mc = monteCarloUni(sampleN,alpha,numArms);

%% Pick usertypes from the data
%     sampleN = samples;
%     sampled_theta = zeros(sampleN, numArms, numTypes) ;
%     sampled_prob = zeros(sampleN,numArms) ;
%     for I = 1 : numArms
%        for J = 1 : numTypes
%            sampled_theta(:,I,J) = random('beta',priors(1,I,J),priors(2,I,J),sampleN,1) ; 
%        end
%     end
%    % Sum s
%         for I = 1 : sampleN
%             sampled_prob(I,:) = sampled_theta(I,:,allUsers(I)) ;
%         end
%     
% %     bestArms = zeros(sampleN,numTypes);
% %     topArm = zeros(1,numTypes);
% %     for I = 1 : sampleN
% %         for J = 1 : numTypes
% %             [~,bestArms(I,J)] = max(sampled_theta(I,:,J));
% %         end
% %     end
% %     final = zeros(numTypes,2);
% %     for I = 1:numTypes
% %         topArm(1,I) = mode(bestArms(:,I),1);
% %         final(I,:) = [topArm(1,I),histc(bestArms(:,I),topArm(1,I))];
% %     end
% 
%     
%     % Calculate top two arms
%     numWins = zeros(numArms) ;
%     [~,armIdx] = max(sampled_prob.') ;
%     for J = 1 : numArms
%        tmp = find(armIdx == J) ;
%        numWins(J,1) = length(tmp) ;
%     end
%     numWins = [(1:numArms).',numWins] ;
%     mc  = numWins(:,2) ./ sampleN;
% %     
% % numWins = sortrows(numWins,-2) ;
% % arm1 = numWins(1,1) ; arm2 = numWins(2,1) ;
% % 
% % twoArmSamples = sampled_prob(:,[arm1 arm2]) ;
% % [~,armIdxTwoOnly] = max(twoArmSamples.') ;
% % tmp = find(armIdxTwoOnly == 1) ;
% % arm1Prob = length(tmp)/sampleN ;
% % [arm1 arm1Prob]

end

