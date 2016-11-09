clear all

load ConversionData ;
% data 7 columns are: time, days, variation, context_1 context 2 
%                   context_3 converted

% plot variations shows anomaly
% plot(convData(:,3),'+') ; % variation 5 is not real
% Remove variation 5 from data
% idx = find(convData(:,3) == 5) ;
% convData(idx,:) = [] ;
% idx = find(convData(:,3) > 5) ;
% convData(idx,3) = convData(idx,3) - 1 ;



% Count number of conversions
totalNumConv = sum(convData(:,7)) ;
ratioConv = totalNumConv/size(convData,1) ;

% Get number of variations
numVariations = length(unique(convData(:,3))) ;

% Histogram the Successes for each Context
% Grab only the successes
idx = find(convData(:,7) == 1) ;
convDataSucc = convData(idx,:) ;
% figure; hist(convDataSucc(:,1))
% figure; hist(convDataSucc(:,2))
% figure; hist(convDataSucc(:,4))
% figure; hist(convDataSucc(:,5))
% figure; hist(convDataSucc(:,6))
% 
% figure;hist3(convDataSucc(:,[1,6]))
% figure;hist3(convDataSucc(:,[2,6]))
% figure;hist3(convDataSucc(:,[4,6]))
% figure;hist3(convDataSucc(:,[5,6]))
% 
% figure; hist3(convDataSucc(:,[4,5]))
% figure; hist3(convDataSucc(:,[4,6]))
% figure; hist3(convDataSucc(:,[5,6]))
% 
% figure; hist3(convDataSucc(:,[2,4]))
% figure; hist3(convDataSucc(:,[2,5]))

%% Develop algorithm to Segment Data into Relevant Classes

% Known Context: Day of Week
% get data
dayData = convData(:,2) ;
hData = histogram(dayData,'BinWidth',1,'BinMethod','integers') ;
dayCounts = hData.Values.' ; dayNum = hData.BinEdges.' + 0.5 ; 
dayNum(end) = [] ; % remove last day because it is at the last edge
close(gcf) ;
% Divide into Week Day and WeekEnd
dayOfWeek = mod(dayNum,7) ;
weekEndIndex = [] ;
for I = 0 : 6
    dayOneIndex = find(dayOfWeek == I) ;
    dayTwoIndex = find(dayOfWeek == mod(I+1,7)) ;
    weekEndIndex{1,I+1} = [dayOneIndex;dayTwoIndex] ;
    weekEndCount(1,I+1) = sum(dayCounts([dayOneIndex;dayTwoIndex])) ;
end
[~,Idx] = min(weekEndCount) ; %Assumption min count occurs on weekend
weekEndDays = sort(dayNum(weekEndIndex{1,Idx}),1,'ascend') ;
% Classify dayData as Week Day or WeekEnd
classData = zeros(size(convData)) ;
for I = 1 : length(convData)
    if (mod(convData(I,2),7) == weekEndDays(1)) || (mod(convData(I,2),7) == weekEndDays(2))
        classData(I,2) = 0 ;
    else
        classData(I,2) = 1 ;
    end
end
    
% Known Context: Time of Day
% get data
timeData = round(convData(:,1)) ;
timeData = mod(timeData,24) ;
hData = histogram(timeData,'BinWidth',1,'BinMethod','integers') ;
timeCounts = hData.Values.' ; timeVal = hData.BinEdges.' + 0.5 ; 
timeVal(end) = [] ; % remove last time because it is at the last edge
close(gcf)
[~,timeIdx] = sort(timeCounts,1,'descend') ;
awakeHrs = timeVal(timeIdx(1:14)) ; % 14 hour awake time
% check if awake hours wraps the 24 hour boundary
wrapCheck = max(awakeHrs) - min(awakeHrs) ;
% Classify timeData as awake or sleep
tmpData = mod(round(convData(:,1)),24) ;
if wrapCheck == 23
    sleepHrs = timeVal(timeIdx(15:24)) ;
    minHr = min(sleepHrs) ; maxHr = max(sleepHrs) ;
    classData(find(tmpData < minHr | tmpData > maxHr),1) = 1 ;
else
    minHr = min(awakeHrs) ; maxHr = max(awakeHrs) ;
    classData(find(tmpData >= minHr & tmpData <= maxHr),1) = 1 ;
end

% Classify Context 1
% get Data
CData = convDataSucc(:,4) ;
hData = histogram(CData,'BinWidth',1,'BinMethod','integers') ;
dataCounts = hData.Values.' ; dataVal = hData.BinEdges.' + 0.5 ; 
dataVal(end) = [] ; % remove last time because it is at the last edge
close(gcf) ;
[~,dataIdx] = sort(dataCounts,1,'descend') ;
bigDataVal = dataVal(dataIdx(1:6)) ; % by eye
% Classify CData as Big or Small
for I = 1 : length(bigDataVal)
    classData( find(convData(:,4) == bigDataVal(I)), 4 ) = 1 ;
end

% Classify Context 2
% get Data
CData = convDataSucc(:,5) ;
hData = histogram(CData,'BinWidth',1,'BinMethod','integers') ;
dataCounts = hData.Values.' ; dataVal = hData.BinEdges.' + 0.5 ; 
dataVal(end) = [] ; % remove last time because it is at the last edge
close(gcf) ;
[~,dataIdx] = sort(dataCounts,1,'descend') ;
bigDataVal = dataVal(dataIdx(1:7)) ; % by eye
% Classify CData as Big or Small
for I = 1 : length(bigDataVal)
    classData( find(convData(:,5) == bigDataVal(I)), 5 ) = 1 ;
end

% Classify Context 3
% get Data
CData = convDataSucc(:,6) ;
hData = histogram(CData,'BinWidth',1,'BinMethod','integers') ;
dataCounts = hData.Values.' ; dataVal = hData.BinEdges.' + 0.5 ; 
dataVal(end) = [] ; % remove last time because it is at the last edge
close(gcf) ;
[~,dataIdx] = sort(dataCounts,1,'descend') ;
bigDataVal = dataVal(dataIdx(1:2)) ; % by eye
% Classify CData as Big or Small
classData( find(convData(:,6) == bigDataVal(1)), 6 ) = 2 ;
classData( find(convData(:,6) == bigDataVal(2)), 6 ) = 1 ;
