function Ks = Kbulk_s__(x,y)

faults=1;

% sG = miu;
% return
Ks =nan(size(x));

%%% Sarah Bunney
%layer A
[ai,bi]=find(y>-600);
Ks(ai,bi)=30e9;
% sG(ai,bi)=4.10e9;
%layer B
[ai,bi]=find(y<=-600 & y>-2700);
Ks(ai,bi)=30e9;
%layer C
[ai,bi]=find(y<=-2700 & y>-4000);
Ks(ai,bi)=30e9; 
%layer D
[ai,bi]=find(y<=-4000 & y>-7500);
Ks(ai,bi)=100e9;
%layer E
[ai,bi]=find(y<=-7500 & y>-9000);
Ks(ai,bi)=30e9;
%layer F
[ai,bi]=find(y<=-9000);
Ks(ai,bi)=100e9;


%%% Sarah Bunney + Alia (densities taken from Alia up to 1.5 Km)
% % % %layer A
% % % [ai,bi]=find(y>-600);
% % % sG(ai,bi)=1.17e9;
% % % %layer B
% % % [ai,bi]=find(y<=-600 & y>-1100);
% % % sG(ai,bi)=6.09e9;
% % % %layer C
% % % [ai,bi]=find(y<=-1100 & y>-1500);
% % % sG(ai,bi)=6.76e9;
% % % %layer D
% % % [ai,bi]=find(y<=-1500 & y>-2700);
% % % sG(ai,bi)=7.79e9;
% % % %layer E
% % % [ai,bi]=find(y<=-2700 & y>-4000);
% % % sG(ai,bi)=16.26e9;
% % % %layer F
% % % [ai,bi]=find(y<=-4000 & y>-7500);
% % % sG(ai,bi)=25.10e9;
% % % %layer G
% % % [ai,bi]=find(y<=-7500 & y>-9000);
% % % sG(ai,bi)=4.87e9;
% % % %layer H
% % % [ai,bi]=find(y<=-9000);
% % % sG(ai,bi)=33.87e9;


% % % if faults
% % % %     %     Fault1P1=[4000 -1500];
% % % %     Fault1P1=[2000 -372.5-3725];
% % % %     Fault1P2=[3000 -372.5];
% % % %     %     Fault2P1=[7600 -1500];
% % % %     Fault2P1=[7260 -3000];
% % % %     Fault2P2=[7940 0];
% % %     
% % %     Fault1P1=[3770.77 -1500];
% % %     Fault1P2=[4000 -200];
% % %     Fault2P1=[7098.08 -1500];
% % %     Fault2P2=[7500    0];
% % %     
% % %     Fault_core_radius=12;
% % %     Fault_damage_radius=50;
% % %     for i=1:length(x(:))
% % %         P=[x(i) y(i)];
% % %         if OnDislocation(P,Fault1P1,Fault1P2) && abs(phi2_disl_(x(i),y(i),Fault1P1,Fault1P2))<Fault_damage_radius
% % %             if abs(phi2_disl_(x(i),y(i),Fault1P1,Fault1P2))<Fault_core_radius
% % %                 sG(i)=3.84e7;
% % %             else
% % %                 %                 tht=(abs(phi2_disl_(x(i),y(i),Fault1P1,Fault1P2)) - Fault_core_radius)/(Fault_damage_radius-Fault_core_radius);
% % %                 %                 sG(i)=tht*sG(i)+(1-tht)*3.84e7;
% % %                 sG(i)=3.84e8;
% % %             end
% % %         elseif OnDislocation(P,Fault2P1,Fault2P2) && abs(phi2_disl_(x(i),y(i),Fault2P1,Fault2P2))<Fault_damage_radius
% % %             if abs(phi2_disl_(x(i),y(i),Fault2P1,Fault2P2))<Fault_core_radius
% % %                 sG(i)=3.84e7;
% % %             else
% % %                 %                             tht=(abs(phi2_disl_(x(i),y(i),Fault2P1,Fault2P2)) - Fault_core_radius)/(Fault_damage_radius-Fault_core_radius);
% % %                 %                             sG(i)=tht*sG(i)+(1-tht)*3.84e7;
% % %                 sG(i)=3.84e8;
% % %             end
% % %         end
% % %     end
% % % end

