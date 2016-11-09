
clear all

%% CMAB 1.3:
% Construct distributions by following data to create our world and use CMAB to test
% Probabilities updated early hour (not every visitor)

%% Notes
% .* [1 2 4 8 16]) ---> Converts a 5 bit binary to its decimal
% contextUser(I,5) ~= 2 ---> The 3rd entry of context 3 we ignore

load ConversionData ; % variable name is convData
% data 7 columns are: 
% time, days, variation, context_1 context 2 context_3 converted

load FieldClassData ; % variable name is classData

%% Metrics of the Data
numVisitors = size(convData,1) ;            % Length if data
numArms = length(unique(convData(:,3))) ;   % Number of arms
numContexts = size(convData,2) - 2 ;        % Contexts used
numTypes = 2^numContexts ;                  % Total number of users possible
armPrior = [1; 1] ;                         % Prior success & failures
noise = 0.0*randn(1);                       % Gaussian noise

% setup arrays
numTrials = zeros(numVisitors,numArms,numTypes) ;   % Num times armXcontext occurred
numSucc = zeros(numVisitors,numArms,numTypes) ;     % Num times armXcontext succeeded
contextUser = classData ;                           % size runN by numContexts
theta = zeros(numArms,numTypes) ;                   % Conversion rates
numTrialWorld = zeros(numArms,numTypes) ;          % Num times armXcontext occurred in data
numSuccWorld = zeros(numArms,numTypes) ;            % Num times armXcontext succeeded in data

% To update arm probabilities
hours = convData(:,1,:) + (24 .* convData(:,2,:));
numHours = ceil(hours(end));
hourSucc = zeros(numHours,numArms,numTypes);
hourTrials = zeros(numHours,numArms,numTypes);

% To update certainity probabilities
days = convData(:,2,:);

%% Construct distributions of World using field data
% Uses success & failures of arms and their contexts to build "Real World" distributions

for I = 1 : numVisitors
    % get userType
    userType = 1 + sum(contextUser(I,:) .* [1 2 4 8 16]) ;
    if contextUser(I,5) ~= 2
        armPlayedWorld = convData(I,3) + 1 ;
        numTrialWorld(armPlayedWorld,userType) = numTrialWorld(armPlayedWorld,userType) + 1 ;
        numSuccWorld(armPlayedWorld,userType) = ...
            numSuccWorld(armPlayedWorld,userType) + convData(I,7) ;
    end
end

% Grab conversion rates for extra
for I = 1 : numArms
    for J = 1 : numTypes
        theta(I,J) = numSuccWorld(I,J)/numTrialWorld(I,J) ;
    end
end
        
% Strangely userType 17 to 28 is never seen !!!!
% e.g., numPlayedWorld(1-8,userType = 17 to 28) = 0
% we can then reduce the number of userTypes        
        
        
%% Instead of taking theta conv rates, use "Real World" distributions
sampled_theta = zeros(numArms,1) ;
progX = TextProgressBar(0) ;
countI = 0 ;
currentHour = 0;
hourIdx = find(hours <= currentHour + 1);


for I = 1 : numVisitors
    progX = TextProgressBar(I/numVisitors,progX) ;
    if contextUser(I,5) ~= 2
        countI = countI + 1 ;
        % convert context to  a number between 1 and 32
        userType = 1 + sum(contextUser(I,:) .* [1 2 4 8 16]) ;
        allUsers(countI,1) = userType ; % keep track of all userTypes
        for J = 1 : numArms
            curPrior(1,1) = armPrior(1) + hourSucc(currentHour + 1,J,userType) ;
            curPrior(2,1) = armPrior(2) + ...
                hourTrials(currentHour + 1,J,userType) - hourSucc(currentHour + 1,J,userType) ;
            % construct sample
            sampled_theta(J,1) = random('beta',curPrior(1),curPrior(2),1) ;
        end
        [~,armIdx] = max(sampled_theta) ;
        playArm = zeros(1,numArms,numTypes) ;
        % added some noise in the randn signal
        playArm(1,armIdx,userType) = ...
            (rand(1) < random('beta',numSuccWorld(armIdx,userType),numTrialWorld(armIdx,userType),1) + noise) ;
        % update number of successes
        numSucc(countI+1,:,:) = numSucc(countI,:,:) + playArm ;
        % update number of successes
        trialVec = zeros(1,numArms,numTypes) ; trialVec(1,armIdx,userType) = 1 ;
        numTrials(countI+1,:,:) = numTrials(countI,:,:) + trialVec ;
    end
    
    if I == hourIdx(end - 1)
        hourSucc(currentHour+1,:,:) = numSucc(countI+1,:,:);
        hourTrials(currentHour+1,:,:) = numTrials(countI+1,:,:);
        if I ~= numVisitors - 1
            currentHour = currentHour + floor((hours(hourIdx(end)+1))) - floor(hours(hourIdx(end)));
        end
        hourIdx = find(hours < currentHour + 1);
    end
    
end
numRuns = countI ; % Use numRuns instead of runN because one context has
% a 2 which we ignore

TextProgressBar('close');
pause

%% Monte Carlo Sampling
% Find certainity probability of best arm being the best

% Number of samples:
sampleN = 10000 ; 

% Collection of the samples across arm and context
sampled_theta = zeros(sampleN, numArms, numTypes) ;

% Amalgamated samples across arms
sampled_prob = zeros(sampleN,numArms) ;

% Priors of all armXcontext at end
finalPrior(1,:,:) = ones(size(numSucc(1,:,:))) + numSucc(numRuns,:,:) ;
finalPrior(2,:,:) = ones(size(numSucc(1,:,:))) + ...
    numTrials(numRuns,:,:) - numSucc(numRuns,:,:) ;

% findProb = monteCarlo(sampleN, finalPrior, numArms, numTypes, allUsers);

% Sample sampleN times from each Arm
progX = TextProgressBar(0) ;
for I = 1 : numArms
    for J = 1 : numTypes
        sampled_theta(:,I,J) = random('beta',finalPrior(1,I,J),finalPrior(2,I,J),sampleN,1) ;
    end
    progX = TextProgressBar(I/numArms,progX) ;
end
TextProgressBar('close');

% CHECK THIS...
% Sums s
for I = 1 : sampleN
    sampled_prob(I,:) = sampled_theta(I,:,allUsers(I)) ;
end

% Find best arm for each usertype
bestArms = zeros(sampleN,numTypes);
topArm = zeros(1,32);
for I = 1 : sampleN
    for J = 1 : numTypes
        [~,bestArms(I,J)] = max(sampled_theta(I,:,J));
    end
end
final = zeros(numTypes,2);
for I = 1:numTypes
    topArm(1,I) = mode(bestArms(:,I),1);
    final(I,:) = [topArm(1,I),histc(bestArms(:,I),topArm(1,I))];
end

% Use to count [1,histc(bestArms(:,1),1)];

%% Calculate top two arms
numWins = zeros(numArms,1) ;
[~,armIdx] = max(sampled_prob.') ;
for J = 1 : numArms
    tmp = find(armIdx == J) ;
    numWins(J,1) = length(tmp) ;
end
numWins = [(1:numArms).',numWins] ;
allProb = numWins(:,2) ./ sampleN;
numWins = sortrows(numWins,-2) ;
arm1 = numWins(1,1) ; arm2 = numWins(2,1) ;

twoArmSamples = sampled_prob(:,[arm1 arm2]) ;
[~,armIdxTwoOnly] = max(twoArmSamples.') ;
tmp = find(armIdxTwoOnly == 1) ;
arm1Prob = length(tmp)/sampleN ;
[arm1 arm1Prob]

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

betaPrior1 = squeeze( armPrior(1) + numSucc(numVisitors+1,:,:) ) ;
betaPrior2 = squeeze( armPrior(2) + numTrials(numVisitors+1,:,:) ...
    - numSucc(numVisitors+1,:,:) ) ;

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