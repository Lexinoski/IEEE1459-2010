clc
clear
close all

k=3;           %Number of ciclos to integration. (k=15 It allow until 4 refresh per second and inter/subharmonics h=0.0667)
f=60;           %Fundamental Frequency
Fs=25000;       %Sampling Frequency. (Fs=25000 It allow until 208th harmonic.)

T=1/f;          %Fundamental Period.

%% Query Waves of tests. Waves Simulation in Time Domain. Three phase system.
    t=(0:k*T*Fs-1)/Fs;
% Values of prototype developed.
%
            Vef=127;
            Ief=22.727272;
            thetaV=pi*0/180;        %Phase angle to voltage.
            thetaI=pi*0/180;        %Phase angle to current.
            
            va=sqrt(2)*Vef*sin(2*pi*60*t+thetaV);
            vb=sqrt(2)*Vef*sin(2*pi*60*t+thetaV-2*pi/3);
            vc=sqrt(2)*Vef*sin(2*pi*60*t+thetaV+2*pi/3);             
            
            ia=sqrt(2)*Ief*sin(2*pi*60*t+thetaI+pi/6);
            ib=sqrt(2)*Ief*sin(2*pi*60*t+thetaI-2*pi/3-pi/6);
            ic=zeros(size(ia));   
%             
%A.1 (IEEE1459-2010) - The effect of the integration interval:            
%            va=sqrt(2)*0.00035*sin(0.0217*2*pi*60*t+pi*-90/180)+sqrt(2)*0.00014*sin(0.0433*2*pi*60*t+pi*-107.3/180)+sqrt(2)*0.56000*sin(0.9580*2*pi*60*t+pi*-97.2/180)+sqrt(2)*70.7100*sin(1.0000*2*pi*60*t+pi*-7.2/180)+sqrt(2)*0.46000*sin(1.0220*2*pi*60*t+pi*-82.9/180)+sqrt(2)*0.35000*sin(1.0430*2*pi*60*t+pi*-104.3/180)+sqrt(2)*5.02000*sin(3.0000*2*pi*60*t+pi*-76.0/180)+sqrt(2)*0.95000*sin(4.0268*2*pi*60*t+pi*176.4/180)+sqrt(2)*3.18000*sin(5.0000*2*pi*60*t+pi*-114.0/180)+sqrt(2)*2.33000*sin(7.0000*2*pi*60*t+pi*-142.0/180)+sqrt(2)*1.13000*sin(9.0000*2*pi*60*t+pi*-165/180);
%            vb=sqrt(2)*127*sin(2*pi*60*t-2*pi/3);
%            vc=sqrt(2)*127*sin(2*pi*60*t+2*pi/3);
% 
%            ia=sqrt(2)*1.48000*sin(0.0217*2*pi*60*t+pi*80.2/180)+sqrt(2)*2.26000*sin(0.0433*2*pi*60*t+pi*-8.4/180)+sqrt(2)*0.92000*sin(0.9570*2*pi*60*t+pi*-173.5/180)+sqrt(2)*2.24000*sin(0.9580*2*pi*60*t+pi*-193.2/180)+sqrt(2)*70.7100*sin(1.0000*2*pi*60*t+pi*-42.4/180)+sqrt(2)*1.75000*sin(1.0220*2*pi*60*t+pi*-178.9/180)+sqrt(2)*0.91000*sin(1.0430*2*pi*60*t+pi*-202.3/180)+sqrt(2)*19.0900*sin(3.0000*2*pi*60*t+pi*18.30/180)+sqrt(2)*5.43000*sin(4.0268*2*pi*60*t+pi*-87.0/180)+sqrt(2)*7.64000*sin(5.0000*2*pi*60*t+pi*-15.80/180)+sqrt(2)*3.68000*sin(7.0000*2*pi*60*t+pi*-43.2/180)+sqrt(2)*1.41000*sin(9.0000*2*pi*60*t+pi*-69/180);
%            ib=zeros(size(ia));
%            ic=zeros(size(ia));

% Annex B (IEEE1459-2010) - Practical sudies and measurements:
%             va=sqrt(2)*100*sin(1*2*pi*60*t+pi*-0/180)+sqrt(2)*8*sin(3*2*pi*60*t+pi*-70/180)+sqrt(2)*15*sin(5*2*pi*60*t+pi*+140/180)+sqrt(2)*5*sin(7*2*pi*60*t+pi*+20/180);%   
%             vb=sqrt(2)*127*sin(2*pi*60*t-2*pi/3);
%             vc=sqrt(2)*127*sin(2*pi*60*t+2*pi/3);
% 
%             ia=sqrt(2)*100*sin(1*2*pi*60*t+pi*-30/180)+sqrt(2)*20*sin(3*2*pi*60*t+pi*-165/180)+sqrt(2)*15*sin(5*2*pi*60*t+pi*+234/180)+sqrt(2)*10*sin(7*2*pi*60*t+pi*+234/180);
%             ib=zeros(size(ia));
%             ic=zeros(size(ia));            

%% Three-phase nonsinusoidal and unbalanced systems (IEEE Std 1459/2010 - Item 3.2.3 pgs 25 a 28)

%LKC - Law of Kirchhoff of Currents
in=-(ia+ib+ic);         % in=0 to three wire sistems

%Fourier Transform to analysis in the frequency domain.
Vaw=fft(va);            % Voltages Fourier Transform
Vbw=fft(vb);
Vcw=fft(vc);
Vabw=fft(va-vb);
Vbcw=fft(vb-vc);
Vcaw=fft(vc-va);
Iaw=fft(ia);            % Current Fourier Transform
Ibw=fft(ib);
Icw=fft(ic);
Inw=fft(in);

Vaw=[Vaw(1)/(k*T*Fs) sqrt(2)*Vaw(2:((k*T*Fs)/2)+1)/(k*T*Fs)];    %Harmonic Phasors (Rectangular), to h=0 (component DC position 1 center of fft, it do not need modulation. (Effectives Values)
Vbw=[Vbw(1)/(k*T*Fs) sqrt(2)*Vbw(2:((k*T*Fs)/2)+1)/(k*T*Fs)];
Vcw=[Vcw(1)/(k*T*Fs) sqrt(2)*Vcw(2:((k*T*Fs)/2)+1)/(k*T*Fs)];
Vabw=[Vabw(1)/(k*T*Fs) sqrt(2)*Vabw(2:((k*T*Fs)/2)+1)/(k*T*Fs)];
Vbcw=[Vbcw(1)/(k*T*Fs) sqrt(2)*Vbcw(2:((k*T*Fs)/2)+1)/(k*T*Fs)];
Vcaw=[Vcaw(1)/(k*T*Fs) sqrt(2)*Vcaw(2:((k*T*Fs)/2)+1)/(k*T*Fs)];
Iaw=[Iaw(1)/(k*T*Fs) sqrt(2)*Iaw(2:((k*T*Fs)/2)+1)/(k*T*Fs)]; 
Ibw=[Ibw(1)/(k*T*Fs) sqrt(2)*Ibw(2:((k*T*Fs)/2)+1)/(k*T*Fs)]; 
Icw=[Icw(1)/(k*T*Fs) sqrt(2)*Icw(2:((k*T*Fs)/2)+1)/(k*T*Fs)]; 
Inw=[Iaw(1)/(k*T*Fs) sqrt(2)*Inw(2:((k*T*Fs)/2)+1)/(k*T*Fs)]; 

Va=sqrt((k*T*Fs/2+1)*mean(abs(Vaw).^2));  %Total RMS Voltages and Currents.
Vb=sqrt((k*T*Fs/2+1)*mean(abs(Vbw).^2));
Vc=sqrt((k*T*Fs/2+1)*mean(abs(Vcw).^2));
Vab=sqrt((k*T*Fs/2+1)*mean(abs(Vabw).^2));
Vbc=sqrt((k*T*Fs/2+1)*mean(abs(Vbcw).^2));
Vca=sqrt((k*T*Fs/2+1)*mean(abs(Vcaw).^2));
Ia=sqrt((k*T*Fs/2+1)*mean(abs(Iaw).^2));
Ib=sqrt((k*T*Fs/2+1)*mean(abs(Ibw).^2));
Ic=sqrt((k*T*Fs/2+1)*mean(abs(Icw).^2));
In=sqrt((k*T*Fs/2+1)*mean(abs(Inw).^2));

Va1=abs(Vaw(k+1));  %Fundamental RMS Voltages and Currents.
Vb1=abs(Vbw(k+1));
Vc1=abs(Vcw(k+1));

Ia1=abs(Iaw(k+1));
Ib1=abs(Ibw(k+1));
Ic1=abs(Icw(k+1));
In1=abs(Inw(k+1));

%% Calculation of the Powers in Frequency Domain.

%Active Power per phase 
Pa=sum(real(Vaw.*conj(Iaw)));      %Active Power (IEEE Std 1459/2010 - Item 3.1.2.3 - pg 9)
Pb=sum(real(Vbw.*conj(Ibw)));      %Fundamental and harmonics Active Powers is the part real of (Vw.*conj(Iw)).
Pc=sum(real(Vcw.*conj(Icw))); 

%Apparent Powers per phase (Magnitude).
Sa=Va*Ia;                          %Apparent Power(IEEE Std 1459/2010 - Item 3.1.2.7 - pg 10)
Sb=Vb*Ib;
Sc=Vc*Ic;

%Nonactive powers per phase
Na=sign(imag(Vaw(k+1)*conj(Iaw(k+1))))*sqrt(Sa^2-Pa^2); %Nonactive power. (IEEE Std 1459/2010 - Item 3.1.2.14 - pg 12)
Nb=sign(imag(Vbw(k+1)*conj(Ibw(k+1))))*sqrt(Sb^2-Pb^2); %Signal of the Fundamental Reactive Power.
Nc=sign(imag(Vcw(k+1)*conj(Icw(k+1))))*sqrt(Sc^2-Pc^2);

%% Calculation of the Power Factor per phase. (IEEE Std 1459/2010 - Item 3.1.1.5 pg 5)

FPa=Pa/Sa;
FPb=Pb/Sb;
FPc=Pc/Sc;

%% Calculation of S and FP three phase unbalanced. (IEEE Std 1459/2010 - Item 3.2.2 pg 15)

P=Pa+Pb+Pc;         %Active Power Total (IEEE Std 1459/2010 - Equation 3.2.2.2 pg 16)
N=Na+Nb+Nc;         %Reactive Power Total (IEEE Std 1459/2010 - Equation 3.2.2.3 pg 17)

SA=Sa+Sb+Sc;        %Arithmetic apparent power (IEEE Std 1459/2010 - Equation 3.2.2.5 pg 18)
SV=sqrt(P^2+N^2);   %Vetorial apparent power (IEEE Std 1459/2010 - Equation 3.2.2.6 pg 19)

FPV=P/SV;           %Vetorial power factor (IEEE Std 1459/2010 - Equation 3.2.2.7 pgs 19)
FPA=P/SA;           %Arithmetic power factor (IEEE Std 1459/2010 - Equation 3.2.2.7 pgs 20)

%% Calculation of Se and FPe (Effective) three phase. (IEEE Std 1459/2010 - Item 3.2.2.8, 3.2.2.9 and 3.2.3.1 pgs 22 to 29)

% Effectives Voltage end Current, Total and Fundamental

Ve=sqrt((3*(Va^2+Vb^2+Vc^2)+Vab^2+Vbc^2+Vca^2)/18);     %Effective Voltage. (IEEE Std 1459/2010 - Item 3.2.2.8 pgs 23)
Ie=sqrt((Ia^2+Ib^2+Ic^2+In^2)/3);                       %Effective Current. (IEEE Std 1459/2010 - Item 3.2.2.8 pgs 22)

Ve1=sqrt((3*(abs(Vaw(k+1))^2+abs(Vbw(k+1))^2+abs(Vcw(k+1))^2)+abs(Vabw(k+1))^2+abs(Vbcw(k+1))^2+abs(Vcaw(k+1))^2)/18);      %Fundamental Effective Voltage (IEEE Std 1459/2010 - Ítem 3.2.3.1 pgs 27)
Ie1=sqrt((abs(Iaw(k+1))^2+abs(Ibw(k+1))^2+abs(Icw(k+1))^2+abs(Inw(k+1))^2)/3);                                              %Fundamental Effective Current (IEEE Std 1459/2010 - Ítem 3.2.3.1 pgs 26)


Se=3*Ve*Ie;     % Effective apparent power. (IEEE Std 1459/2010 - Equation 3.2.2.8 pg 24)
FPe=P/Se;       % Effective power fector. (IEEE Std 1459/2010 - Equation 3.2.2.9 pg 29)

%% Calculation of the all THDs.

THDva=100*(sqrt(Va^2-Va1^2))/Va1;          %Calculation of the THD per phase (IEEE Std 1459/2010 - Equation 3.1.2.1 pg 8)
THDvb=100*(sqrt(Vb^2-Vb1^2))/Vb1;  
THDvc=100*(sqrt(Vc^2-Vc1^2))/Vc1;  
THDia=100*(sqrt(Ia^2-Ia1^2))/Ia1; 
THDib=100*(sqrt(Ib^2-Ib1^2))/Ib1; 
THDic=100*(sqrt(Ic^2-Ic1^2))/Ic1; 

THDv =100*(sqrt(Ve^2-Ve1^2))/Ve1;          %Calculation of the THDs effectives to three phase system(IEEE Std 1459/2010 - Equation 3.2.3.1 pg 28)
THDi =100*(sqrt(Ie^2-Ie1^2))/Ie1;

%% Show values
fprintf('Pa = %.2f W\n',Pa);
fprintf('Na = %.2f var\n\n',Na);
%fprintf('Va = %.2f V\n',Va);
%fprintf('Ia = %.2f A\n',Ia);
%fprintf('THDva = %.4f \n',THDva);
%fprintf('THDia = %.4f \n\n',THDia);

fprintf('Pb = %.2f W\n',Pb);
fprintf('Nb = %.2f var\n\n',Nb);
 
fprintf('Pc = %.2f W\n',Pc);
fprintf('Nc = %.2f var\n\n',Nc);

fprintf('FPa = %.4f \n',FPa);
fprintf('FPb = %.4f \n',FPb);
fprintf('FPc = %.4f \n\n',FPc);

fprintf('FPV = %.4f \n',FPV);
fprintf('FPA = %.4f \n',FPA);
fprintf('FPe = %.4f \n\n',FPe);

%% Graphic of the Va and Ia in Time Domain.
figure
plot(t,va,'k--','LineWidth',2');
grid on
hold on
plot(t,ia,'r:','LineWidth',2');
legend('Tensão','Corrente');

%% Graphic of the Va and Ia components in Frequency Domain.
figure
fw=stem(t(1:((k*T*Fs)/2)+1)*Fs/(k*T), abs(Vaw));
grid on
xlabel('Frequency (Hz)')
ylabel('Magnitude (V)')
title('Voltage in Frequency Domain')

figure
fw=stem(t(1:((k*T*Fs)/2)+1)*Fs/(k*T), abs(Iaw));
grid on
xlabel('Frequency (Hz)')
ylabel('Magnitude (A)')
title('Current in Frequency Domain')

