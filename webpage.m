
xstep = 100;

length_target = 0.1; %target length
length_source = 0.1; %source length

epsilon0 = 200;
Csource = 10;
Ctarget = 70;
T1source = 1;
T1target = 10;

length_total = length_target+length_source;
Dsource = 5e-4*(Csource/70)^(1/3);
Dtarget = 5e-4*(Ctarget/70)^(1/3);
alpha = 1/sqrt(Dtarget*T1target);
beta = 1/sqrt(Dsource*T1source);

xtarget = linspace(0,length_target,xstep);
xsource = linspace(length_target,length_total,xstep);
x = [xtarget, xsource];

lambda = beta*Csource*Dsource/(alpha*Ctarget*Dtarget);

E = epsilon0-1;
F = lambda*coth(alpha*length_target)*(sinh(beta*length_target)-cosh(beta*length_target)*tanh(beta*length_total));
G = -cosh(beta*length_target)+sinh(beta*length_target)*tanh(beta*length_total);
CIII = E/(F+G);

A = epsilon0-1;
B = CIII*(cosh(beta*length_target)-sinh(beta*length_target)*tanh(beta*length_total));
C = cosh(alpha*length_target);
CI = (A+B)/C;

Ptarget = 1+CI*cosh(alpha.*xtarget);
Etarget = trapz(Ptarget)/xstep

Psource = epsilon0+CIII*(cosh(beta.*xsource)-sinh(beta.*xsource)*tanh(beta*length_total));
Esource = trapz(Psource)/xstep

P = [Ptarget,Psource];

clf(figure(1));figure(1);
plot(x,P,'LineWidth',1.5);
set(gca,'FontSize',22,'LineWidth',1.5);
xlabel('Position (\mum)')
ylabel('Polarization \muwaves on')
grid on;
y1=get(gca,'ylim');
hold on
grid on
plot([length_target length_target],y1,'k','LineWidth',1.5);
hold off  


