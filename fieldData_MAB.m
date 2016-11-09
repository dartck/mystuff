
clear all

load ConversionData ; % variable name is convData
% data 7 columns are: 
% time, days, variation, context_1 context 2 context_3 converted

load FieldClassData ; % variable name is classData

% Type of simulation
% simType = 0 ; % No Context MAB, Using calculated fixed conversion rates
%               for each arm
% simType = 1 ; % Contextual MAB, Using calculated fixed conversion 
%                 probabiities with no noise for each context
% simType = 2 ; % No Context MAB, Let the real data choose the arm !
simType = 3 ; % Contextual MAB, Let the real data choose the arm !

% Metrics of the Data
runN = size(convData,1) ;
numArms = length(unique(convData(:,3))) ; % 8 in this case
numContexts = size(convData,2) - 2 ; % 5 in this case
numTypes = 2^numContexts ; % two classes per context
armPrior = [1; 1] ;

switch simType
    case {0} % no context MAB
        numTrials = zeros(runN,numArms) ;
        numSucc = zeros(runN,numArms) ;
        contextUser = classData ;
        
        % Extract conversion rates from Field data
        theta = zeros(numArms,1) ;
        numPlayed = zeros(numArms,1) ;
        numSuccTrials = zeros(numArms,1) ;
        for I = 1 : runN
            % get userType
            userType = 1 + sum(contextUser(I,:) .* [1 2 4 8 16]) ;
            if contextUser(I,5) ~= 2
                armPlayed = convData(I,3) + 1 ;
                numPlayed(armPlayed,1) = numPlayed(armPlayed,1) + 1 ;
                numSuccTrials(armPlayed,1) = ...
                    numSuccTrials(armPlayed,1) + convData(I,7) ;
            end
        end
        for I = 1 : numArms
            theta(I,1) = numSuccTrials(I,1)/numPlayed(I,1) ;
        end
        
    case {1} % contextual MAB with fixed conversion rates
        % setup arrays
        numTrials = zeros(runN,numArms,numTypes) ;
        numSucc = zeros(runN,numArms,numTypes) ;
        % contextArm = zeros(numArms,numContexts) ;
        contextUser = classData ; % size runN by numContexts
        
        
        % Extract conversion rates from Field data
        theta = zeros(numArms,numTypes) ;
        numPlayed = zeros(numArms,numTypes) ;
        numSuccTrials = zeros(numArms,numTypes) ;
        for I = 1 : runN
            % get userType
            userType = 1 + sum(contextUser(I,:) .* [1 2 4 8 16]) ;
            if contextUser(I,5) ~= 2
                armPlayed = convData(I,3) + 1 ;
                numPlayed(armPlayed,userType) = numPlayed(armPlayed,userType) + 1 ;
                numSuccTrials(armPlayed,userType) = ...
                    numSuccTrials(armPlayed,userType) + convData(I,7) ;
            end
        end
        for I = 1 : numArms
            for J = 1 : numTypes
                theta(I,J) = numSuccTrials(I,J)/numPlayed(I,J) ;
            end
        end
        % Strangely userType 17 to 28 is never seen !!!!
        % e.g., numPlayed(1-8,userType = 17 to 28) = 0
        % we can then reduce the number of userTypes
        
    case {2} % no context MAB, let the real data choose the Arm
        numTrials = zeros(runN,numArms) ;
        numSucc = zeros(runN,numArms) ;
        contextUser = classData ;
        
        % Extract conversion rates from Field data
        theta = zeros(numArms,1) ;
        numPlayed = zeros(numArms,1) ;
        numSuccTrials = zeros(numArms,1) ;
        for I = 1 : runN
            % get userType
            userType = 1 + sum(contextUser(I,:) .* [1 2 4 8 16]) ;
            if contextUser(I,5) ~= 2
                armPlayed = convData(I,3) + 1 ;
                numPlayed(armPlayed,1) = numPlayed(armPlayed,1) + 1 ;
                numSuccTrials(armPlayed,1) = ...
                    numSuccTrials(armPlayed,1) + convData(I,7) ;
            end
        end
        for I = 1 : numArms
            theta(I,1) = numSuccTrials(I,1)/numPlayed(I,1) ;
        end

    case {3} % contextual MAB letting real data choose Arm
        % setup arrays
        numTrials = zeros(runN,numArms,numTypes) ;
        numSucc = zeros(runN,numArms,numTypes) ;
        % contextArm = zeros(numArms,numContexts) ;
        contextUser = classData ; % size runN by numContexts
        
        
        % Extract conversion rates from Field data
        theta = zeros(numArms,numTypes) ;
        numPlayed = zeros(numArms,numTypes) ;
        numSuccTrials = zeros(numArms,numTypes) ;
        for I = 1 : runN
            % get userType
            userType = 1 + sum(contextUser(I,:) .* [1 2 4 8 16]) ;
            if contextUser(I,5) ~= 2
                armPlayed = convData(I,3) + 1 ;
                numPlayed(armPlayed,userType) = numPlayed(armPlayed,userType) + 1 ;
                numSuccTrials(armPlayed,userType) = ...
                    numSuccTrials(armPlayed,userType) + convData(I,7) ;
            end
        end
        for I = 1 : numArms
            for J = 1 : numTypes
                theta(I,J) = numSuccTrials(I,J)/numPlayed(I,J) ;
            end
        end
        % Strangely userType 17 to 28 is never seen !!!!
        % e.g., numPlayed(1-8,userType = 17 to 28) = 0
        % we can then reduce the number of userTypes        
        
end


% Setup Conversion rates
% convScaleFactor = 0.1 ;
% theta = zeros(numArms,numTypes) ;
% for I = 1 : 2^numContexts
%         mystr = dec2bin(I-1) ;
%         lenStr = length(mystr) ;
%         contextVec = zeros(1,numContexts) ;
%         for K = 1 : lenStr
%             contextVec(K) = str2num(mystr(K)) ;
%         end
%         for J = 1 : numArms
%             tmpVal = convScaleFactor*sum((contextArm(J,:) .* contextVec))/numContexts ;
%             mnVal = tmpVal*ones(60,1) ;
%             sigma = 0.1 * ones(60,1) ;
%             theta(J,I) =  max(0,mean(normrnd(mnVal,sigma))) ;
%         end
%         fullContext(I,:) = contextVec ;
% end
% if sum(sum(theta > 1)) > 1
%     error('Error in Theta\n')
% end

switch simType
    case{0}
        % Playing the arms
        sampled_theta = zeros(numArms,1) ;
        % Run the bandit
        progX = TextProgressBar(0) ;
        countI = 0 ;
        for I = 1 : runN
            progX = TextProgressBar(I/runN,progX) ;
            if contextUser(I,5) ~= 2
                countI = countI + 1 ;
                % convert context to  a number between 1 and 32
                userType = 1 + sum(contextUser(I,:) .* [1 2 4 8 16]) ;
                allUsers(countI,1) = userType ; % keep track of all userTypes
                for J = 1 : numArms
                    curPrior(1,1) = armPrior(1) + numSucc(countI,J) ;
                    curPrior(2,1) = armPrior(2) + ...
                        numTrials(countI,J) - numSucc(countI,J) ;
                    % construct sample
                    sampled_theta(J,1) = random('beta',curPrior(1),curPrior(2),1) ;
                end
                [~,armIdx] = max(sampled_theta) ;
                playArm = zeros(1,numArms) ;
                % added some noise in the randn signal
                playArm(1,armIdx) = (rand(1) < theta(armIdx) ) ;
                % update number of successes
                numSucc(countI+1,:) = numSucc(countI,:) + playArm ;
                % update number of successes
                trialVec = zeros(1,numArms) ; trialVec(1,armIdx) = 1 ;
                numTrials(countI+1,:) = numTrials(countI,:) + trialVec ;
            end
            
        end
        numRuns = countI ; % Use numRuns instead of runN because one context has
        % a 2 which we ignore
        
        TextProgressBar('close');
        
    case{1}
        % Playing the arms
        sampled_theta = zeros(numArms,1) ;
        % Run the bandit
        progX = TextProgressBar(0) ;
        countI = 0 ;
        for I = 1 : runN
            progX = TextProgressBar(I/runN,progX) ;
            if contextUser(I,5) ~= 2
                countI = countI + 1 ;
                % convert context to  a number between 1 and 32
                userType = 1 + sum(contextUser(I,:) .* [1 2 4 8 16]) ;
                allUsers(countI,1) = userType ; % keep track of all userTypes
                for J = 1 : numArms
                    curPrior(1,1) = armPrior(1) + numSucc(countI,J,userType) ;
                    curPrior(2,1) = armPrior(2) + ...
                        numTrials(countI,J,userType) - numSucc(countI,J,userType) ;
                    % construct sample
                    sampled_theta(J,1) = random('beta',curPrior(1),curPrior(2),1) ;
                end
                [~,armIdx] = max(sampled_theta) ;
                playArm = zeros(1,numArms,numTypes) ;
                % added some noise in the randn signal
                playArm(1,armIdx,userType) = (rand(1) < theta(armIdx,userType) ) ;
                % update number of successes
                numSucc(countI+1,:,:) = numSucc(countI,:,:) + playArm ;
                % update number of successes
                trialVec = zeros(1,numArms,numTypes) ; trialVec(1,armIdx,userType) = 1 ;
                numTrials(countI+1,:,:) = numTrials(countI,:,:) + trialVec ;
            end
            
        end
        numRuns = countI ; % Use numRuns instead of runN because one context has
        % a 2 which we ignore
        
        TextProgressBar('close');
        
    case{2}
        % Playing the arms
        % Run the bandit
        progX = TextProgressBar(0) ;
        countI = 0 ;
        for I = 1 : runN
            progX = TextProgressBar(I/runN,progX) ;
            if contextUser(I,5) ~= 2
                countI = countI + 1 ;
                % convert context to  a number between 1 and 32
                userType = 1 + sum(contextUser(I,:) .* [1 2 4 8 16]) ;
                allUsers(countI,1) = userType ; % keep track of all userTypes
                playArm = zeros(1,numArms) ;
                playArm(1,convData(countI,3)+1) = convData(countI,7) ;
                % update number of successes
                numSucc(countI+1,:) = numSucc(countI,:) + playArm ;
                % update number of successes
                trialVec = zeros(1,numArms) ; trialVec(1,convData(countI,3)+1) = 1 ;
                numTrials(countI+1,:) = numTrials(countI,:) + trialVec ;
            end
            
        end
        numRuns = countI ; % Use numRuns instead of runN because one context has
        % a 2 which we ignore
        
        TextProgressBar('close');

    case{3}
        % Playing the arms
        % Run the bandit
        progX = TextProgressBar(0) ;
        countI = 0 ;
        for I = 1 : runN
            progX = TextProgressBar(I/runN,progX) ;
            if contextUser(I,5) ~= 2
                countI = countI + 1 ;
                % convert context to  a number between 1 and 32
                userType = 1 + sum(contextUser(I,:) .* [1 2 4 8 16]) ;
                allUsers(countI,1) = userType ; % keep track of all userTypes
                playArm = zeros(1,numArms,numTypes) ;
                playArm(1,convData(countI,3)+1,userType) = convData(countI,7) ;
                % update number of successes
                numSucc(countI+1,:,:) = numSucc(countI,:,:) + playArm ;
                % update number of trials
                trialVec = zeros(1,numArms,numTypes) ; trialVec(1,convData(countI,3)+1,userType) = 1 ;
                numTrials(countI+1,:,:) = numTrials(countI,:,:) + trialVec ;
            end
            
        end
        numRuns = countI ; % Use numRuns instead of runN because one context has
        % a 2 which we ignore
        
        TextProgressBar('close');
        
end


%%%% Analysis of Simulation %%%%%%%%%%%%%

%% SimType0
% plot trial results
% figure ;
% plot(numTrials(1:numRuns,:)) ;
% % conversion rates
% figure ;
% convRate = numSucc(1:numRuns,:) ./ numTrials(1:numRuns,:) ;
% plot(convRate(:,:))
% theta

%% SimType1
% plot trial results for specific user type
% figure ;
% userType = 10 ;
% plot(numTrials(1:numRuns,:,userType)) ;
% 
% % plot conversion rates
% figure ;
% userType  = 4 ;
% convRate = numSucc(1:numRuns,:,:) ./ numTrials(1:numRuns,:,:) ;
% plot(convRate(:,:,userType))




% questionable at 1/16
sampleN = numRuns ;

switch simType
    case {0}
        pause
    case{1}
        % calculate arm providing best conversion rate for the given population
        for I = 1 : numTypes
            numUsers(I,1) = length(find(allUsers == I)) ;
        end
        for I = 1 : numRuns
            convRate = squeeze( numSucc(I,:,:) ./ numTrials(I,:,:) ) ;
            armConvNum = sum(bsxfun(@times,convRate,numUsers.'),2) ;
            [~,bestArmIdx(I,1)] = max(armConvNum) ;
        end
        
        for I = 1 : 32
            userType = I ;
            [I;theta(:,I)]
            % plot trial results for specific user type
            figure ;
            plot(numTrials(1:numRuns,:,userType)) ;
            
            % plot conversion rates
            figure ;
            convRate = numSucc(1:numRuns,:,:) ./ numTrials(1:numRuns,:,:) ;
            plot(convRate(:,:,userType))
            pause
            close all
        end
    case{2}
        % calculae probability that arm A is better than arm B
        % First sample sampleN times from each beta distribution
        
        sampled_theta = zeros(sampleN,numArms) ;
        
        % Calculate beta distribution priors at the end
        for J = 1 : numArms
            curPrior(1,J) = armPrior(1) + numSucc(numRuns,J) ;
            curPrior(2,J) = armPrior(2) + ...
                numTrials(numRuns,J) - numSucc(numRuns,J) ;
        end
        
        % Sample sampleN times from each Arm
        for J = 1 : numArms
            sampled_theta(:,J) = random('beta',curPrior(1,J),curPrior(2,J),sampleN,1) ;
        end
        
        % Calculate top two
        [~,armIdx] = max(sampled_theta.') ;
        for J = 1 : numArms
            tmp = find(armIdx == J) ;
            numWins(J,1) = length(tmp) ;
        end
        numWins = [(1:numArms).',numWins] ;
        numWins = sortrows(numWins,-2) ;
        arm1 = numWins(1,1) ; arm2 = numWins(2,1) ;
        
        twoArmSamples = sampled_theta(:,[arm1 arm2]) ;
        [~,armIdx] = max(twoArmSamples.') ;
        tmp = find(armIdx == 1) ;
        arm1Prob = length(tmp)/sampleN ;
        [arm1 arm1Prob]

    case{3}
        % calculae probability that arm A is better than arm B
        % First sample sampleN times from each beta distribution
        
        % Calculate beta distribution priors at the end
        
        sampled_theta = zeros(sampleN, numArms, numTypes) ;
        
        curPrior(1,:,:) = ones(size(numSucc(1,:,:))) + numSucc(numRuns,:,:) ;
        curPrior(2,:,:) = ones(size(numSucc(1,:,:))) + ...
                          numTrials(numRuns,:,:) - numSucc(numRuns,:,:) ;
        
        % Sample sampleN times from each Arm
        progX = TextProgressBar(0) ;
        for I = 1 : numArms
            for J = 1 : numTypes
                sampled_theta(:,I,J) = random('beta',curPrior(1,I,J),curPrior(2,I,J),sampleN,1) ;
            end
            progX = TextProgressBar(I/numArms,progX) ;
        end
        TextProgressBar('close');
        
        sampled_prob = zeros(sampleN,numArms) ;
        for I = 1 : sampleN
            sampled_prob(I,:) = sampled_theta(I,:,allUsers(I)) ;
        end
        
        
        % Calculate top two
        numWins = zeros(numArms,1) ;
        [~,armIdx] = max(sampled_prob.') ;
        for J = 1 : numArms
            tmp = find(armIdx == J) ;
            numWins(J,1) = length(tmp) ;
        end
        numWins = [(1:numArms).',numWins] ;
        numWins = sortrows(numWins,-2) ;
        arm1 = numWins(1,1) ; arm2 = numWins(2,1) ;
        
        twoArmSamples = sampled_prob(:,[arm1 arm2]) ;
        [~,armIdxTwoOnly] = max(twoArmSamples.') ;
        tmp = find(armIdxTwoOnly == 1) ;
        arm1Prob = length(tmp)/sampleN ;
        [arm1 arm1Prob]
        
%%%%%%% All Wrong because every user Type sampled equally %%%%%%        
%         % Calculate top two
%         % max across userTypes
%         [maxVals, typeIdx] = max(sampled_theta,[],3) ; % 3rd col is userType
%         % max across the arms for each trial
%         [~,armWin] = max(maxVals.') ;
%         armWin = armWin.' ; % armIdx shows which arm won
%         % which user Type won
%         for I = 1 : sampleN
%             typeWin(I,1) = typeIdx(I,armWin(I)) ;
%         end
%  
%         for J = 1 : numArms
%             tmp = find(armWin == J) ;
%             numWins(J,1) = length(tmp) ;
%         end
%         numWins = [(1:numArms).',numWins] ;
%         numWins = sortrows(numWins,-2) ;
%         arm1 = numWins(1,1) ; arm2 = numWins(2,1) ;
%         
%         twoArmSamples = sampled_theta(:,[arm1 arm2],:) ;
%         [maxVals, twoTypeIdx] = max(twoArmSamples,[],3) ; % 3rd col is userType
%         [~,twoArmWin] = max(maxVals.') ;
%         twoArmWin = twoArmWin.' ; % armIdx shows which arm won       
%         tmp = find(armWin == 1) ;
%         arm1Prob = length(tmp)/sampleN ;
%         [arm1 arm1Prob]        
        
end

pause
%%%% We have examined to here %%%%%%%%%

%% Linear Regression to predict Arm contexts

% predition model theta_hat(n) = aCoefVec (dot product) userContextVec(n) +
%                                 phi
for J = 1 : numArms
    
    for I = 1 : numContexts
        aCoefVec(I,1:numContexts) = (fullContext(:,I).') * fullContext ;
        phiCoefVec(I,1) = sum(fullContext(:,I)) ;
        constVec(I,1) = (fullContext(:,I).') * theta(J,:).' ;
    end
    
    aCoefVec(numContexts+1,1:numContexts) = sum(fullContext,1) ;
    phiCoefVec(numContexts+1,1) = numTypes ;
    constVec(numContexts+1,1) = sum(theta(J,:)) ;
    
    coefMat = [phiCoefVec aCoefVec] ;
    armContext(J,:) = (coefMat \ constVec).' ;
    
    armContextScale(J,1) = mean(contextArm(J,:) ./ armContext(J,2:numContexts+1)) ;
    
end

% Scale the predicted arm context (for comparison)
scaledArmContext = bsxfun(@times,armContext,armContextScale) ;
for J = 1 : numArms
    figure ;
    plot(contextArm(J,:),scaledArmContext(J,2:numContexts+1),'+') ;
end

%% Calculating the probabilities

betaPrior1 = squeeze( armPrior(1) + numSucc(runN+1,:,:) ) ;
betaPrior2 = squeeze( armPrior(2) + numTrials(runN+1,:,:) ...
    - numSucc(runN+1,:,:) ) ;

Npoints = 128 ;
xVal = linspace(0,1,Npoints).' ;
permuteIdx = [1 2 3;2 1 3;3 2 1] ;
probDist = zeros(Npoints,numArms) ;

for J = 1 : numTypes
    %%%%%%% we hard-code three arm probability sphere here %%%%%%%%%%%%%
    probSphere = ones(Npoints,Npoints,Npoints) ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for I = 1 : numArms
        tmpProb = [] ;
        probDist(:,I) = betapdf(xVal, betaPrior1(I,J), betaPrior2(I,J) ) ;
        probDist(:,I) = probDist(:,I)/sum(probDist(:,I)) ;
        tmpProb(1:Npoints,:,:) = probDist(:,I) ;
        tmpProb = permute(tmpProb,permuteIdx(I,:)) ;
        probSphere = probSphere .* ...
            bsxfun(@times,ones(Npoints,Npoints,Npoints),tmpProb) ;
    end
    
    % Calculate probability each Arm is better
    for I = 1 : numArms
        tmpProbSphere = permute(probSphere,permuteIdx(I,:)) ;
        tmpArmProb = 0 ; tmpExpectCR = 0 ;
        for K = 1 : Npoints
            tmpArmProb = tmpArmProb + ...
                sum(sum(sum(tmpProbSphere(K,1:K,1:K)))) ;
            tmpExpectCR = tmpExpectCR + ...
                sum(sum(sum(tmpProbSphere(K,:,:) * xVal(K)))) ;
        end
        armProb(I,J) = tmpArmProb ;
        expectCR(I,J) = tmpExpectCR ;
        expectNum(I,J) = tmpExpectCR * numUsers(J) ;
    end
    armProb(:,J) = armProb(:,J)/sum(armProb(:,J)) ;
end
armProb
expectNum
%% Calculate Arm Probabilities across Contexts
% As number of cases increases exponentially with number of Contexts
% keep number of cases less than a reasonable amount and recursively
% work through the data

armProbPrev = armProb ;
expectNumPrev = expectNum ;
numTypesRecur = numTypes ;
maxNumCases = 1e4 ;
numConGroups = floor(log(maxNumCases)/log(numArms)) ;
numCases = numArms^numConGroups ;
while numTypesRecur > 1

    totalArmProb = [] ;
    totalExpNum = [] ;
    numSubDiv = floor(numTypesRecur/numConGroups) ;
    
    if numSubDiv > 0
        probVal = ones(numCases,numSubDiv) ;
        numConvPerArm = zeros(numCases,numArms,numSubDiv) ;
        for I = 1 : numSubDiv
            for J = 1 : numCases
                mystr = dec2base(J-1,numArms) ;
                lenStr = length(mystr) ;
                caseVec = zeros(1,numConGroups) ;
                for K = 1 : lenStr
                    caseVec(K) = str2num(mystr(K)) ;
                end
                for K = 1 : numConGroups
                    probVal(J,I) = probVal(J,I) * ...
                        armProbPrev(caseVec(K)+1,(I-1)*numConGroups+K) ;
                    numConvPerArm(J,caseVec(K)+1,I) = ...
                        numConvPerArm(J,caseVec(K)+1,I) + ...
                        expectNumPrev(caseVec(K)+1,(I-1)*numConGroups+K) ;
                end
            end
            probVal(:,I) = probVal(:,I)/sum(probVal(:,I)) ;
        end
        
        % Add up Probabilities across SubDivisions
        totalArmProb = zeros(numArms,numSubDiv) ;
        totalExpNum = zeros(numArms,numSubDiv) ;
        for I = 1 : size(numConvPerArm,3)
            [~,maxIdx] = max(numConvPerArm(:,:,I),[],2) ;
            for J = 1 : numArms
                tmpIdx = find(maxIdx == J) ;
                totalArmProb(J,I) = sum(probVal(tmpIdx,I)) ;
                totalExpNum(J,I) = sum(probVal(tmpIdx,I) .* numConvPerArm(tmpIdx,J,I)) ;
            end
        end
        
    end
    
    
    if rem(numTypesRecur,numConGroups) > 0
    
        numConGroupsLast = numTypesRecur - numSubDiv*numConGroups ;
        numCasesLast = numArms^numConGroupsLast ;
        numConvPerArmLast = zeros(numCasesLast,numArms) ;
        probValLast = ones(numCasesLast, 1) ;
        for J = 1 : numCasesLast
            mystr = dec2base(J-1,numArms) ;
            lenStr = length(mystr) ;
            caseVec = zeros(1,numConGroupsLast) ;
            for K = 1 : lenStr
                caseVec(K) = str2num(mystr(K)) ;
            end
            for K = 1 : numConGroupsLast
                probValLast(J,1) = probValLast(J,1) * ...
                    armProbPrev(caseVec(K)+1,numSubDiv*numConGroups+K) ;
                numConvPerArmLast(J,caseVec(K)+1) = ...
                    numConvPerArmLast(J,caseVec(K)+1) + ...
                    expectNumPrev(caseVec(K)+1,(I-1)*numConGroups+K) ;
            end
        end
        probValLast(:,1) = probValLast(:,1)/sum(probValLast) ;
        [~,maxIdx] = max(numConvPerArmLast,[],2) ;
        for J = 1 : numArms
            tmpIdx = find(maxIdx == J) ;
            armProbLast(J,1) = sum(probValLast(tmpIdx)) ;
            expNumLast(J,1) = sum(probValLast(tmpIdx) .* numConvPerArmLast(tmpIdx,J)) ;
        end
        totalArmProb(:,numSubDiv+1) = armProbLast ;
        totalExpNum(:,numSubDiv+1) = expNumLast ;
    end
    armProbPrev = totalArmProb ;
    expNumPrev = totalExpNum ;
    numTypesRecur = size(armProbPrev,2) ;
    
end
format long
finalArmProb = armProbPrev 
finalExpNum = expNumPrev 