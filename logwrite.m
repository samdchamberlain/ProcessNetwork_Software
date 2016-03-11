function [msg] = logwrite(msg,vout)
% Appends a message to the processLog, which is a cell of strings
% Note: to access processLog in the main workspace, type >> global processLog
%
% ------ Inputs ------
% msg = a character array with the message to log to processLog
% vout = output message to screen as well? (0 = no, 1 = yes)
%
% ------ Outputs -----
% msg = the message logged in processLog
% 
% --------------------

global processLog

if ~iscell(processLog)
    processLog = {};     
end

if ~ischar(msg)
    warning('Warning: Processing message to log is not a string. Message not written.')
else
    processLog{size(processLog,1)+1,1} = msg;
end
    
if vout
    disp(msg)
end
    