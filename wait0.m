
for iter=1:1000
    load coupling_with_MUFITS/SCENARIO12_with_matlab//WAIT.txt
    while WAIT==0
        pause(1)
        load WAIT.txt
    end
    %do something
    WAIT=0;
    fid = fopen('WAIT.txt', 'wt');
    fprintf(fid,'%d',WAIT);
    fclose(fid);
end