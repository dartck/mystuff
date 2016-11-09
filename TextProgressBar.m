function newVal = TextProgressBar(newVal,prvVal)
% TEXTPROGRESSBAR Text-based progress bar.
%
% Normal usage:
%
%     x = TextProgressBar(0) ;  % initialisation
%     for I = 0:.01:1 
%         x = TextProgressBar(I,x) ;  % update
%         pause(.1) ; 
%     end 
%     TextProgressBar('close') ;  % close the bar
%
% Note that the input values (progress) are expected to be between 0 and 1.
%
% The bar can also be used backwards, from 1 to 0:
%     x = TextProgressBar(1) ;  % initialisation at 100% 
%     for I = 0:.01:1 
%         x = TextProgressBar(1-I,x) ;  % update
%         pause(.1) ; 
%     end 
%     TextProgressBar('close') ;  % close the bar


% Action = 1 -> initialise the progress bar
% Action = 2 -> update the bar
% Action = 3 -> close the bar
if nargin < 2
    prvVal = 0 ;
    if ischar(newVal)
        action = 3 ;
    else
        action = 1 ;
    end
else
    action = 2 ;    
end

% Make sure the input values are between 0 and 1
prvVal = min(max(prvVal,0),1) ;
newVal = min(max(newVal,0),1) ;

% Round "percentage" values
prvPrc = round(prvVal*50) ;
newPrc = round(newVal*50) ;

% Initialise/Update/Close the progress bar
switch action
    case 1
        % Initialise the progress bar
        fprintf( ...
            ['   [' repmat('%%',1,newPrc) ...
            repmat('_',1,50-newPrc) '] ' ...
            repmat(' ',1,3-length(num2str(2*newPrc))) ...
            num2str(2*newPrc) '%%']) ;
    case 2
        % If the percent value is the same, do nothing
        if newPrc == prvPrc
            newVal = newPrc/50 ;
            return
            % Else, update
        elseif newPrc > prvPrc
            fprintf([ ...
                repmat('\b',1,6+(50-prvPrc)) ...
                repmat('%%',1,(newPrc-prvPrc)) ...
                repmat('_',1,50-newPrc) '] ' ...
                repmat(' ',1,3-length(num2str(2*newPrc))) ...
                num2str(2*newPrc) '%%']) ;
        else
            fprintf([ ...
                repmat('\b',1,6+(50-newPrc)) ...
                repmat('_',1,50-newPrc) '] ' ...
                repmat(' ',1,3-length(num2str(2*newPrc))) ...
                num2str(2*newPrc) '%%']) ;
        end
    case 3
        % Close the Progress bar (start a new line)
        fprintf('\n') ;
end

% Output: "new value" in percent
newVal = newPrc/50 ;


