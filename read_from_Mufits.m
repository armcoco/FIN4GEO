function [startRow, endRow] = read_from_Mufits(filename)


%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.

dataArray{1}='SOMETHING_DIFFERENT_FROM_DATA';
startRow=1;
while ~strcmp(dataArray{1},'DATA')
    dataArray=textscan(fileID,'%s%[^\n\r]',1);
    startRow=startRow+1;
end    
endRow=startRow;
while ~strcmp(dataArray{1},'ENDDATA')
    dataArray=textscan(fileID,'%s%[^\n\r]',1);
    endRow=endRow+1;
end    
endRow=endRow-3;
fclose(fileID);

