function dataArray=read_data_from_Mufits(filename,startRow,endRow,ncols)

% formatSpec = '%f%f%f%[^\n\r]';
formatSpec='';
for i=1:ncols
formatSpec = [formatSpec, '%f'];
end
formatSpec = [formatSpec, '%[^\n\r]'];

fileID = fopen(filename,'r');
textscan(fileID, '%[^\n\r]', startRow-1, 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.