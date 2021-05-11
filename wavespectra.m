%%WAVE SPECTRA plot ---- spectral density versus frequency ratio
%%Jonswap spectrum
hs=5; %wave height in m
t=10; %time period in seconds
v=3; %peakedness factor chosen between 1 to 7
g=9.81; %gravitational constant
w=0:0.1:3; %frequency is the varying component
n=length(w);
wo=(2*pi)/t;
alpha=3.25*(10^-3)*(hs^2)*(wo^4)*(1-(0.287*(log(v))));
for i=1:n
if w(i)<=wo
sigma(i)=0.07;%spectral width parameter
else
sigma(i)=0.09;
end
x(i)=-((w(i)-wo)^2)/(2*(sigma(i)^2)*((wo)^2));
y(i)=-1.25*((w(i)/wo)^(-4));
aw(i)=exp(x(i));
z(i)=exp(y(i))*(v^aw(i))*alpha*(g^2);
s(i)=z(i)/(w(i)^5);
p(i)=w(i)/wo;
i=i+1;
end
%%PM spectrum
hs=5; %wave height in m
t=10; %time period in seconds
g=9.81; %gravitational constant
v=20; %mean wind speed in m/s
wo=(2*pi)/t;
w=0:0.1:3; %frequency is the varying component
n=length(w);
for i=1:n
x(i)=-1.25*((w(i)/wo)^(-4));
a(i)=exp(x(i));
b(i)=1/((w(i))^5);
s1(i)=0.0081*a(i)*b(i)*(g)^2;
p1(i)=w(i)/wo;
i=i+1;
end
%%Modified PM spectrum
hs=5; %wave height in m
t=10; %time period in seconds
wo=(2*pi)/t;
w=0:0.1:3; %frequency is the varying component
n=length(w);
for i=1:n
y(i)=(-1.25)*((w(i)/wo)^(-4));
a(i)=exp(y(i));
b(i)=(wo^4)/((w(i))^5);
s2(i)=0.3125*((hs^2)*a(i)*b(i));
p2(i)=w(i)/wo;
i=i+1;
end
%%ISSC spectrum
hs=5; %wave height in m
t=10; %time period in seconds
wo=(2*pi)/t;
w=0:0.1:3; %frequency is the varying component
n=length(w);
for i=1:n
y(i)=-1.2489*((w(i)/wo)^(-4));
a(i)=exp(y(i));
b(i)=(wo^4)/((w(i))^5);
s3(i)=0.3123*((hs^2)*a(i)*b(i));
p3(i)=w(i)/wo;
i=i+1;
end
plot(p,s,'k'); %jonswap spectrum
hold on;
plot(p1,s1,'r');%Modified PM spectrum
hold on;
plot(p2,s2,'b'); %Bretschneider spectrum
hold on;
plot(p3,s3,'g'); %ISSC spectrum
xlabel('Frequency Ratio');
ylabel('Spectral density');
title('WAVE SPECTRA (Mean wind speed=20m/s, Wave Height=5m, Time period=10s)')
len1 = length(p1)
interval = (s1(len1)-s1(1))/len1 
bracket = sum(p1)+sum(p1)-p1(1)-p1(len1)
int = interval*bracket/2