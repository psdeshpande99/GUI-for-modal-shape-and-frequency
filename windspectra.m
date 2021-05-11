%%wind spectra plot---- spectral density versus theta
%%davenport spectrum
um=20; %mean wind speed at a height of 10 m
del=0.001; %surface drag coefficient
lu=1200; %integral length for davenport spectrum in m
w=0.001:0.001:10; %frequency is the varying component
theta=(w*lu)/(2*pi*um);
a=4*((theta).^2);
b=(1+(theta.^2)).^(4/3);
x=a./b;
y=(x*del*(um^2));
su=y./w;
%%harris spectrum
um=20; %mean wind speed at aheight of 10 m
del=0.001; %surface drag coefficient
lu=1800; %integral length for harris spectrum in m
w=0.001:0.001:10; %frequency is the varying component
theta1=(w*lu)/(2*pi*um);
a=4*theta1;
b=(2+(theta1.^2)).^(5/6);
x=a./b;
y=x*del*(um^2);
su1=y./w;
%%Kaimal spectrum
t=12; %time period in seconds
w=0.01:0.001:100; %frequency is the varying component
wp=(2*pi)/t; %frequency in radians per second
z=10; %reference height is 10 m
zs=20; %the surface height usually taken as 20 m
uz=(wp*z)/0.025;
if z<=zs
sigma=0.15*uz*((zs/z)^0.125);
else
sigma=0.15*uz*((zs/z)^0.275);
end
theta2=w./wp;
a=6.8*theta2;
b=(1+(10.2*theta2)).^(5/3);
x=a./b;
y=x.*(sigma^2);
su2=y./w;
%%API(2000) spectrum
t=12; %time period in seconds
w=0.01:0.001:100; %frequency is the varying component
wp=(2*pi)/t; %frequency in radians per second
z=10; %reference height is 10 m
zs=20; %the surface height usually taken as 20 m
uz=(wp*z)/0.025;
if z<=zs
sigma=0.15*uz*((zs/z)^0.125);
else
sigma=0.15*uz*((zs/z)^0.275);
end
theta3=w./wp;
b=1+(1.5*(theta3)).^(5/3);
x=theta3./b;
y=x.*(sigma^2);
su3=y./w;
loglog(theta,su,'r','linewidth',2);%Davenport Spectrum
hold on;
loglog(theta1,su1,'b','linewidth',2);%Harris spectrum
hold on;
loglog(theta2,su2,'k', 'linewidth',2);%Kaimal
hold on;
loglog(theta3,su3,'g','linewidth',2);
xlabel('Derivable variable');
ylabel('Spectral density');
title('WIND SPECTRA (Mean wind speed=20m/s, Reference Height=10m, Time period=12s)');
legend('Davenport Spectra','Harris Spectra','Kaimal Spectra','API(2000)Spectra');
