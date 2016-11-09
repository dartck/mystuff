 
clear all

%% CMAB 1.311:
% Simulator version
% Run regular CMAB and extract World thetas
% Follows the data
% But tests arms against 1st CMAB's thetas

load ConversionData ; % variable name is convData
% data 7 columns are: 
% time, days, variation, context_1 context 2 context_3 converted

load FieldClassData ; % variable name is classData
load FieldClassDataSim ; % variable name is classDataSim


%% Metrics of the Data
numVisitors = size(convData,1) ;            % Length if data
numArms = length(unique(convData(:,3))) ;   % Number of arms
numContexts = size(convData,2) - 2 ;        % Contexts used
numTypes = 2^numContexts ;                  % Total number of users possible
armPrior = [1; 1] ;                         % Prior success & failures
noise = 0.0*randn(1);                       % Gaussian noise
sampleN = 1000;                            % Number of samples for Monte Carlo
cutoff = 36;                                % Days of useful data

% setup arrays
numTrialsSim = zeros(numVisitors,numArms,numTypes) ;   % Num times armXcontext occurred in 1st CMAB
numSuccSim = zeros(numVisitors,numArms,numTypes) ;     % Num times armXcontext succeeded in 1st CMAB
numTrials = zeros(numVisitors,numArms,numTypes) ;   % Num times armXcontext occurred
numSucc = zeros(numVisitors,numArms,numTypes) ;     % Num times armXcontext succeeded
contextUser = classData;                           % Data for CMAB 1 in binary
contextUserSim = classDataSim;                      % Data for CMAB 2 in binary
theta = zeros(numArms,numTypes) ;                   % Conversion rates
numTrialWorld = zeros(numArms,numTypes) ;          % Num times armXcontext occurred in data
numSuccWorld = zeros(numArms,numTypes) ;            % Num times armXcontext succeeded in data

% To update arm probabilities per hour
hours = convData(:,1,:) + (24 .* convData(:,2,:));
numHours = ceil(hours(end));
hourSucc = zeros(numHours,numArms,numTypes);
hourTrials = zeros(numHours,numArms,numTypes);

%% Create world distributions for 1st CMAB
% Follow data and extract # plays and # wins for each arm across each
% context
progX = TextProgressBar(0) ;
for I = 1 : numVisitors
    progX = TextProgressBar(I/numVisitors,progX) ;
    % get userType
    userType = 1 + sum(contextUser(I,:) .* [1 2 4 8 16]) ;
    if contextUser(I,5) ~= 2
        armPlayedWorld = convData(I,3) + 1 ;
        numTrialWorld(armPlayedWorld,userType) = numTrialWorld(armPlayedWorld,userType) + 1 ;
        numSuccWorld(armPlayedWorld,userType) = ...
            numSuccWorld(armPlayedWorld,userType) + convData(I,7) ;
    end
end
TextProgressBar('close');
% Finding conversions rates for each arm only (for MAB testing)
armConversions  = sum(numSuccWorld,2);
armPlays = sum(numTrialWorld,2);
MABtheta = armConversions ./ armPlays;

%% Run 1st CMAB and extract thetas for 2nd CMAB

display('Running first CMAB...');
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
        numSuccSim(countI+1,:,:) = numSuccSim(countI,:,:) + playArm ;
        % update number of successes
        trialVec = zeros(1,numArms,numTypes) ; trialVec(1,armIdx,userType) = 1 ;
        numTrialsSim(countI+1,:,:) = numTrialsSim(countI,:,:) + trialVec ;
    end
    
    if I == hourIdx(end - 1)
        hourSucc(currentHour+1,:,:) = numSuccSim(countI+1,:,:);
        hourTrials(currentHour+1,:,:) = numTrialsSim(countI+1,:,:);
        if I ~= numVisitors - 1
            currentHour = currentHour + floor((hours(hourIdx(end)+1))) - floor(hours(hourIdx(end)));
        end
        hourIdx = find(hours < currentHour + 1);
    end   
end

TextProgressBar('close');
numRuns = countI ;
% Extract thetas
armPlaysCont = squeeze(numTrialsSim(numRuns,:,:));
armSuccCont = squeeze(numSuccSim(numRuns,:,:));
CMABtheta = armSuccCont ./ armPlaysCont ;

%% Run 2nd CMAB using data and previous CMAB thetas
% Reset parameters
display('Running second CMAB...');
sampled_theta = zeros(numArms,1) ;
progX = TextProgressBar(0) ;
countI = 0 ;
currentHour = 0;
hourIdx = find(hours <= currentHour + 1);
hourSucc = zeros(numHours,numArms,numTypes);
hourTrials = zeros(numHours,numArms,numTypes);
findProbCont = zeros(numHours,numArms);

for I = 1 : numVisitors
    progX = TextProgressBar(I/numVisitors,progX) ;
    if contextUserSim(I,5) ~= 2
        countI = countI + 1 ;
        % convert context to  a number between 1 and 32
        userType = 1 + sum(contextUserSim(I,:) .* [1 2 4 8 16]) ;
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
            (rand(1) < CMABtheta(armIdx,userType) + noise) ;
        % update number of successes
        numSucc(countI+1,:,:) = numSucc(countI,:,:) + playArm ;
        % update number of successes
        trialVec = zeros(1,numArms,numTypes) ; trialVec(1,armIdx,userType) = 1 ;
        numTrials(countI+1,:,:) = numTrials(countI,:,:) + trialVec ;
    end
    
    if I == hourIdx(end - 1)
        hourSucc(currentHour+1,:,:) = numSucc(countI+1,:,:);
        hourTrials(currentHour+1,:,:) = numTrials(countI+1,:,:);
        findProbCont(currentHour + 1,:) = monteCarloCont(sampleN,...
            [1 + hourSucc(currentHour+1,:,:);1+ hourTrials(currentHour+1,:,:) - hourSucc(currentHour+1,:,:)], numArms, numTypes, allUsers);
        
        if I ~= numVisitors - 1
            currentHour = currentHour + floor((hours(hourIdx(end)+1))) - floor(hours(hourIdx(end)));
        end
        hourIdx = find(hours < currentHour + 1);
    end
    
end
numRuns = countI ; % Use numRuns instead of runN because one context has
% a 2 which we ignore
convRateCont = hourSucc ./hourTrials;

TextProgressBar('close');
figure;
plot(findProbCont);

%pause

%% Regular MAB
%%% Multi-Armed Bandit %%%
display('Running regular MAB...');
progX = TextProgressBar(0) ;
% setup arrays
numTrialsUni = zeros(numRuns,numArms) ;
numSuccUni = zeros(numRuns,numArms) ;
findProbUni = zeros(numHours,numArms);

currentHour = 0;
hourIdx = find(hours <= currentHour + 1);
hourSuccUni = zeros(numHours,numArms);
hourTrialsUni = zeros(numHours,numArms);
% Use conversion rates from contextual bandits

% Playing the arm
sampled_thetaUni = zeros(numArms,1) ;

% Run the bandit
for I = 1 : numVisitors
    progX = TextProgressBar(I/numVisitors,progX) ;
    for J = 1 : numArms
        curPriorUni(1,1) = armPrior(1) + hourSuccUni(currentHour + 1,J) ;
        curPriorUni(2,1) = armPrior(2) + ...
                   hourTrialsUni(currentHour + 1,J) - hourSuccUni(currentHour + 1,J) ;
        % construct sample
        sampled_thetaUni(J,1) = random('beta',curPriorUni(1),curPriorUni(2),1) ;
    end
    [~,armIdxUni] = max(sampled_thetaUni) ;
    playArmUni = zeros(1,numArms) ;
    
    % added some noise in the randn signal
    playArmUni(1,armIdxUni) = (rand(1) < ( MABtheta(armIdxUni) + noise )) ;
    % regretUni(I+1) = regretUni(I) + bestArm(userType) - MABtheta(armIdxUni) ;
    % update number of successes
    numSuccUni(I+1,:) = numSuccUni(I,:) + playArmUni ;
    % update number of trials
    trialVecUni = zeros(1,numArms) ; trialVecUni(1,armIdxUni) = 1 ;
    numTrialsUni(I+1,:) = numTrialsUni(I,:) + trialVecUni ; 
    
    if I == hourIdx(end - 1)
        hourSuccUni(currentHour+1,:) = numSuccUni(I,:);
        hourTrialsUni(currentHour+1,:) = numTrialsUni(I,:);
        findProbUni(currentHour+1,:) = monteCarloUni(sampleN, 1 + ...
            [1 + numSuccUni(currentHour+1,:); 1 + numTrialsUni(currentHour+1,:) - numSuccUni(currentHour+1,:)], numArms);        
        if I ~= numVisitors - 1
            currentHour = currentHour + floor((hours(hourIdx(end)+1))) - floor(hours(hourIdx(end)));
        end
        hourIdx = find(hours < currentHour + 1);
    end
    
end

convRateUni = hourSuccUni./hourTrialsUni;

TextProgressBar('close');

% plot(findProbUni);

pause
%% Monte Carlo Sampling
% Find certainity probability of best arm being the best

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
%% Notes
% How to look for best arm...
% 1) Rank the usertypes for each arm starting from best conversion
% 2) Find each usertype's fraction of traffic
% 3) Weight the contribution of each usertype depending on its amount of traffic
% 4) 
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


% change final
% final = final(:,1:2) ;
% for I=1:32
% bindata{I,1} = dec2bin(I-1,5) ;
% end
% for I=1:32
% textstr = bindata{I,1} ;
% newbindata{I,1} = textstr([1,2,3,4,5]) ;
% end
% newdec = bin2dec(newbindata) + 1 ;
% final(:,3) = newdec ;
% newfinal = sortrows(final,3) ;

%% Plot results
% probUniPlot = [(1:size(findProbUni,1)).'/size(findProbUni,1),findProbUni];
% probContPlot = [(1:size(findProbCont,1)).'/size(findProbCont,1),findProbCont];
% 
figure;
plot(findProbCont);

figure;
plot(findProbUni);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          