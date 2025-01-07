FileSaturation='D:\ThermoPoroElastic\Hydrotherm\HalfSpace\Out_Saturation'

fid = fopen(FileSaturation, 'r');

stringa=fscanf(fid,'%s %s %s %s/n')
R = fscanf(fid, '%g', [10 2*5])
stringa=fscanf(fid,'%s %s %s %s/n')
Z = fscanf(fid, '%g', [10 2*3])

r=[];
for i=2:2:10
    r=[r;R(:,i)];
end

z=[];
for i=2:2:6
    z=[z;Z(:,i)];
end

stringa1=fscanf(fid,'%s %s %s %s %s %s %s %s /n')
stringa2=fscanf(fid,'%s %s %s %s %s %s %s /n')
a=length(strfind(stringa1,'Completed'))

Sat=zeros()
while ~a
for i=1:5
R1 = fscanf(fid, '%g', [1 10 ])
Val=fscanf(fid, '%g', [11 30 ])
Sat()
end

stringa1=fscanf(fid,'%s %s %s %s %s %s %s %s /n')
stringa2=fscanf(fid,'%s %s %s %s %s %s %s /n')
a=length(strfind(stringa1,'Completed'))

end

fclose(fid);



%Density
FileDensity='D:\ThermoPoroElastic\Hydrotherm\HalfSpace\Out_Density'

fid = fopen(FileDensity, 'r');

stringa=fscanf(fid,'%s %s %s %s/n')
R = fscanf(fid, '%g', [10 2*5])
stringa=fscanf(fid,'%s %s %s %s/n')
Z = fscanf(fid, '%g', [10 2*3])

stringa1=fscanf(fid,'%s %s %s %s %s %s %s %s /n')
stringa2=fscanf(fid,'%s %s %s %s %s %s %s /n')
a=length(strfind(stringa1,'Completed'))

while ~a
    for i=1:5
        R1 = fscanf(fid, '%g', [1 10 ])
        Val1=fscanf(fid, '%g', [11 30 ])
    end
    
    stringa=fscanf(fid,'%s %s %s %s %s %s/n')
    
    for i=1:5
        R1 = fscanf(fid, '%g', [1 10 ])
        Val2=fscanf(fid, '%g', [11 30 ])
    end
    
    stringa1=fscanf(fid,'%s %s %s %s %s %s %s %s /n')
    stringa2=fscanf(fid,'%s %s %s %s %s %s %s /n')
    a=length(strfind(stringa1,'Completed'))
    
end

fclose(fid);