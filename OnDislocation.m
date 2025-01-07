function val=OnDislocation(P,P1,P2)

val=(P-P1)*(P2-P1)'>=0 & (P-P2)*(P1-P2)'>=0;

% val=1;