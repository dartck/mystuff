function mc = monteCarloUni( samples, priors, numArms)
%MONTECARLO Summary of this function goes here
%   Detailed explanation goes here

    % Number of samples:
    sampleN = samples ;

    % Collection of the samples across arm and context
    sampled_theta = zeros(sampleN, numArms) ;

    for I = 1 : numArms
       sampled_theta(:,I) = random('beta',priors(1,I),priors(2,I),sampleN,1) ;
    end
    
%     % Sum s
%     for I = 1 : sampleN
%         sampled_prob(I,:) = sampled_theta(I,:,allUsers(I)) ;
%     end
    
    % Calculate probability of best arm 
    numWins = zeros(numArms,1) ;
    [~,armIdx] = max(sampled_theta.') ;
    for J = 1 : numArms
       tmp = find(armIdx == J) ;
       numWins(J,1) = length(tmp) ;
    end
    numWins = [(1:numArms).',numWins] ;
    mc  = numWins(:,2) ./ sampleN;
%     
% numWins = sortrows(numWins,-2) ;
% arm1 = numWins(1,1) ; arm2 = numWins(2,1) ;
% 
% twoArmSamples = sampled_prob(:,[arm1 arm2]) ;
% [~,armIdxTwoOnly] = max(twoArmSamples.') ;
% tmp = find(armIdxTwoOnly == 1) ;
% arm1Prob = length(tmp)/sampleN ;
% [arm1 arm1Prob]

end

