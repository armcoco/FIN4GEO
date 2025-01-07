uchamber=uchamber(:)';
vchamber=vchamber(:)';
figure
hold on
plot(xt__([Ghost.Phi2(:).xc]),yt__([Ghost.Phi2(:).yc]),'*b')
plot(xt__([Ghost.Phi2(:).xc])+uchamber,yt__([Ghost.Phi2(:).yc])+vchamber,'*r')

for ii=1:length(Ghost.Phi2)
    plot([xt__([Ghost.Phi2(ii).xc]),xt__([Ghost.Phi2(ii).xc])+uchamber(ii)],[yt__([Ghost.Phi2(ii).yc]),yt__([Ghost.Phi2(ii).yc])+vchamber(ii)],'-k')
end