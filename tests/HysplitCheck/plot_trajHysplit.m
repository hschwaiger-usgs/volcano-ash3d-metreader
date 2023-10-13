clear all;

tj1=load('ftraj1.dat'); tj1x=tj1(:,1); tj1y=tj1(:,2);
tj2=load('ftraj2.dat'); tj2x=tj2(:,1); tj2y=tj2(:,2);
tj3=load('ftraj3.dat'); tj3x=tj3(:,1); tj3y=tj3(:,2);
tj4=load('ftraj4.dat'); tj4x=tj4(:,1); tj4y=tj4(:,2);
tj5=load('ftraj5.dat'); tj5x=tj5(:,1); tj5y=tj5(:,2);
tj6=load('ftraj6.dat'); tj6x=tj6(:,1); tj6y=tj6(:,2);
tj7=load('ftraj7.dat'); tj7x=tj7(:,1); tj7y=tj7(:,2);

hy1=load('htraj1.dat'); hy1x=hy1(:,1); hy1y=hy1(:,2);
hy2=load('htraj2.dat'); hy2x=hy2(:,1); hy2y=hy2(:,2);
hy3=load('htraj3.dat'); hy3x=hy3(:,1); hy3y=hy3(:,2);
hy4=load('htraj4.dat'); hy4x=hy4(:,1); hy4y=hy4(:,2);
hy5=load('htraj5.dat'); hy5x=hy5(:,1); hy5y=hy5(:,2);
hy6=load('htraj6.dat'); hy6x=hy6(:,1); hy6y=hy6(:,2);
hy7=load('htraj7.dat'); hy7x=hy7(:,1); hy7y=hy7(:,2);

% Make sure Hysplit and traj plot in the same domain
if (hy1x(1)<0.0)
 hy1x = hy1x + 360.0;
 hy2x = hy2x + 360.0;
 hy3x = hy3x + 360.0;
 hy4x = hy4x + 360.0;
 hy5x = hy5x + 360.0;
 hy6x = hy6x + 360.0;
 hy7x = hy7x + 360.0;
end


figure;
hold on;
plot(hy1x,hy1y,'ro-');
plot(hy2x,hy2y,'go-');
plot(hy3x,hy3y,'bo-');
plot(hy4x,hy4y,'mo-');
plot(hy5x,hy5y,'co-');
plot(hy6x,hy6y,'yo-');
plot(hy7x,hy7y,'ko-');

plot(tj1x,tj1y,'r*--');
plot(tj2x,tj2y,'g*--');
plot(tj3x,tj3y,'b*--');
plot(tj4x,tj4y,'m*--');
plot(tj5x,tj5y,'c*--');
plot(tj6x,tj6y,'y*--');
plot(tj7x,tj7y,'k*--');

legend('Hy 5k','Hy 10k','Hy 15k','Hy 20k','Hy 30k','Hy 40k','Hy 50k','Trj 5k','Trj 10k','Trj 15k','Trj 20k','Trj 30k','Trj 40k','Trj 50k')
xlabel('Longitude')
ylabel('Latitude')
title('Volc2')

print -djpg volc1.jpg

hold off;

