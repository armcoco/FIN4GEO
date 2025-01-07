function tt=loadMufitsTimes(filename)

Nfiles=0;
filenum=0;
filename1 = [filename,'.0000.SUM'];
while exist(filename1, 'file') == 2
    Nfiles = Nfiles+1;
    filenum=filenum+1;
    str_filenum=num2str(filenum,'%.4d');
    filename1 = [filename,'.',str_filenum,'.SUM'];
end

tt=zeros(Nfiles,1);
for filenum=(1:Nfiles)-1
    str_filenum=num2str(filenum,'%.4d');
    filename1 = [filename,'.',str_filenum,'.SUM'];
    fileID = fopen(filename1,'r');
    dataArray{1}='SOMETHING_DIFFERENT_FROM_TIME';
    while ~strcmp(dataArray{1},'TIME')
        dataArray=textscan(fileID,'%s%[^\n\r]',1);
    end
    timeArray=textscan(fileID,'%.16f',1);
    tt(filenum+1)= timeArray{1};
    fclose(fileID);
end
