alpha=linspace(0,4500,1000);
z=linspace(-6500.3,-6500,1000);
dalpha=alpha(2)-alpha(1);
dz=z(2)-z(1);
drho=6.2567e+03;
GG=6.67408e-11;

[Alpha,Z]=meshgrid(alpha,z);
dGi=-2*pi*GG.*drho.*dalpha.*dz.*Alpha.*Z./sqrt(Alpha.^2+Z.^2).^3;
dGfinal=sum(dGi(:))*1e8;
disp(['dG is: ',num2str(dGfinal)])
