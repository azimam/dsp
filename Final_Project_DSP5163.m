%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Final Project DSP/5163
% Optimal Wiener Filtration
% Instructor: Dr. Grigoryan
% Azima Motaghi
% Spring 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% Reading  512 sample of the signal from the file "boli.sig"
clc
close all
clear all
fid=fopen('boli.sig','rb'); 
O=fread(fid,'float');
fclose(fid); 
clear fid;
O=O'; 
f=O(1:512);
clear O
N=length(f);
%% 
% 1-Process the signal f(n) degraded by the following smooth filter
%  h = [· · · , 0, 1, 2, 3, 3, 2, 1, 1, · · ·]/13     
% changing the filter arrangment to have (h(0) = 3)

h=[3/13 3/13 2/13 1/13 1/13 zeros(1,N-8) 0 1/13 2/13];
F = fft(f);
H = fft(h);
G = F.*H;
g = ifft(G);
% Calculating RMSE of Smoothed Signal vs. Original Signal
err0 = sqrt((f-g)*(f-g)')/N;
fprintf('\n error1 is %6.4f  \n',err0)
% ploting the smoothed and original signal
figure(1)
subplot(2,1,1);
plot(f,'r');
axis([0 512 -1 40]);
title('the Original Signal');
subplot(2,1,2);
plot(g,'g');
axis([0 512 -1 40]);
title('Smooth filtering, h[n]=[ 0, 1, 2, 3, 3, 2, 1, 1]/13');
%text(200,35,['RMSE = ', num2str(err0)],'Color',[0 0 1]);

%%
% 2-Generating normal distribution noise-signals, nk(n), k = 1 : 4, 
%   with mean 0 and standard deviation .02
nk=zeros(4,N);
for i=1:4
%n=0.02*randn(1,N);
n=normrnd(0,2,1,N);
nk(i,:)=n;
end
 
%%
% 3. Adding the noise-signals to the blur signal, to compose four noisy signals
%    gk(n) = g(n) + nk(n) = (f * h)(n) + nk(n), n = 0 : (N-1), k = 1 : 4
gk=zeros(4,N);
cstring='kbrcmygrm'; % color string
for i=1:4
gk(i,:)=nk(i,:)+g;
figure(2)
subplot(4,1,i)
plot(gk(i,:),cstring(i))
axis([0 512 -10 40]);
title([num2str(i) 'th ' 'Degraded signal with smoth filter  plus noise' ])
figure(3)
subplot(4,1,i)
plot(nk(i,:),cstring(i))
axis([0 512 -10 10]);
title([num2str(i) 'th ' 'Noise Signal' ])
end


%%
% 4a. Designing an inverse filter YInv 
Yinv = conj(H)./(H.*conj(H)+0.0001);


%%
%  4b. Applying  the inverse filter YInv to the noisy signals gk(n)
%**************************************************************************
GK=zeros(4,N);
gk_inv=zeros(4,N);

for i=1:4
GK(i,:)=fft(gk(i,:));
GK_inv(i,:) = Yinv.*GK(i,:);
gk_inv(i,:)=ifft(GK_inv(i,:));
figure(4+i)
subplot(5,1,1)
plot(f,cstring(i))
axis([0 512 -1 40]);
title( 'The Original Signal' )
subplot(5,1,2)
plot(g,cstring(i+1))
axis([0 512 -1 40]);
title( 'Smooth filtering' )
subplot(5,1,3)
plot(nk(i,:),cstring(i+2))
axis([0 512 -10 10]);
title([num2str(i) 'th ' 'Noise ' ])
subplot(5,1,4)
plot(gk_inv(i,:),cstring(i+3))
axis([0 512 -50 100]);
title([num2str(i) 'th ' 'Degraded signal with smoth filter  plus noise' ])
subplot(5,1,5)
plot(gk_inv(i,:),cstring(i+4))
axis([0 512 -60 100]);
title([num2str(i) 'th ' 'Inverse Filtered signal' ])
%  Calculating RMSE of Inverse filter Signal vs. Original Signal   
err_Inv = sqrt((gk_inv(i,:)-f)*(gk_inv(i,:)-f)')/N;
fprintf('\n RMSE of Inverse filter Signal is %f  ',err_Inv)

end



%%
% 5a. Design one Wiener
%**************************************************************************
% Energy Spectrum of Noise & Signal & SNR
%**************************************************************************

for i=1:4
NK (i,:)= fft(nk(i,:));
P_NK(i,:) = NK(i,:).*conj(NK(i,:));
end

P_G = G.*conj(G)
ExpN=(P_NK(1,:)+P_NK(2,:)+P_NK(3,:)+P_NK(4,:) )/4;
SNR =ExpN./(P_G)
for i=1:4
 
Ywin = conj(H)./(conj(H).*H+ 1./SNR(1,:));

% 4b. Apply the inverse filter YInv to the noisy signals gk(n) and calculate
% the error of estimation of the signal f(n)
gk_win(i,:)=ifft(Ywin.*GK(i,:));

figure(8+i)
subplot(5,1,1)
plot(f,cstring(i))
axis([0 512 -1 40]);
title( 'The Original Signal of lentgh ' )

subplot(5,1,2)
plot(g,cstring(i+1))
axis([0 512 -1 40]);
title( 'Smooth filtering' )

subplot(5,1,3)
plot(nk(i,:),cstring(i+2))
axis([0 512 -10 10]);
title([num2str(i) 'th ' 'Noise ' ])

subplot(5,1,4)
plot(gk_inv(i,:),cstring(i+3))
axis([0 512 -1 40]);
title([num2str(i) 'th ' 'inverse filtered signal' ])

subplot(5,1,5)
plot(gk_win(i,:),cstring(i+4))
axis([0 512 -1 40]);
title([num2str(i) 'th ' 'Wiener filtration' ])

%  Calculating RMSE of Wiener filter Signal vs. Original Signal    
err_win = sqrt((gk_win(i,:)-f)*(gk_win(i,:)-f)')/N;
fprintf('\n RMSE of Wiener vs. Original Signal is %f',err_win)
end





%%
for i=1:4
Yhom=sqrt(Yinv.*Ywin);
gk_hom(i,:)=ifft(Yhom.*GK(i,:));


figure(8+i)
subplot(5,1,1)
plot(f,cstring(i))
axis([0 512 -1 40]);
title( 'The Original Signal' )

subplot(5,1,2)
plot(g,cstring(i+1))
axis([0 512 -1 40]);
title( 'Smooth filtering' )

subplot(5,1,3)
plot(nk(i,:),cstring(i+2))
axis([0 512 -10 10]);
title([num2str(i) 'th ' 'Noise ' ])

subplot(5,1,4)
plot(gk(i,:),cstring(i+3))
axis([0 512 -10 40]);
title([num2str(i) 'th ' 'Degraded signal with smoth filter plus noise' ])

subplot(5,1,5)
plot(gk_hom(i,:),cstring(i+4))
axis([0 512 -1 40]);
title([num2str(i) 'th ' 'Homo Filtered signal' ])

%  Calculating RMSE of homomorphic filter Signal
err_hom = sqrt((gk_hom(i,:)-f)*(gk_hom(i,:)-f)')/N;
fprintf('\n RMSE of homomorphic filter Signal is %f  ',err_hom)
end
%%

for i=1:4

figure(12+i)
subplot(6,1,1)
plot(f,cstring(i))
axis([0 512 -1 40]);
title( 'The Original Signal' )

subplot(6,1,2)
plot(g,cstring(i+1))
axis([0 512 -1 40]);
title( 'Smooth filtering' )

subplot(6,1,3)
plot(gk(i,:),cstring(i+2))
axis([0 512 -10 40]);
title([num2str(i) 'th ' 'Degraded signal with smoth filter plus noise ' ])

subplot(6,1,4)
plot(gk_inv(i,:),cstring(i+3))
axis([0 512 -60 100]);
title([num2str(i) 'th ' 'inversed filtered signal' ])

subplot(6,1,5)
plot(gk_win(i,:),cstring(i+4))
axis([0 512 -10 50]);
title([num2str(i) 'th ' 'Wiener filtration' ])

subplot(6,1,6)
plot(gk_hom(i,:),cstring(i+5))
axis([0 512 -10 50]);
title([num2str(i) 'th ' 'Homo Filtered signal' ])
end
%%
%*************************************************************************
% Amplitude Correction Function
%*************************************************************************

   
figure(30+i)
alpha=(H.*conj(H))./(H.*conj(H)+1./SNR);
plot(alpha)
axis tight;
s_title=sprintf('Amplitude Correction Function, \\alpha(w)');
h_title=title(s_title);
set(h_title,'Color',[0 0 1]);


%%

h_new = zeros(1,512);
h_new(254) = 1/13;h_new(255) = 2/13;h_new(256) = 3/13; h_new(257) = 3/13; h_new(258) = 2/13; h_new(259) = 1/13; h_new(260) = 1/13;
H_new = fft(h_new);

figure(35+i)
HdB = -20*log2(abs(H_new))/N;
subplot(5,2,[1 2]);
plot(HdB);
axis([1 512 0 .3]);
title('Absolute value of the transfer function H(w)');
ylabel('-20*log(H)/N','Color',[1 0 1]);


Phase_H = phase(H_new);
subplot(5,1,2);
plot(Phase_H,'g');
axis([1 512 -5 5]);
title('Phase of Transfer Function H(w)');
ylabel('radians','Color',[1 0 1]);

Log_Yinv = -20*log2(abs(Yinv))/N;
subplot(5,1,3);
plot(Log_Yinv,'r');
axis([1 512 -0.2 .3]);
title('Absolute Value of Transfer Function Y_i_n_v(w)');
ylabel('-20*log(H)/N','Color',[1 0 1]);
subplot(5,1,4);
plot(phase(Yinv),'m');
axis([1 512 -5 5]);
title('Phase of Transfer Function Y_i_n_v(w)');
ylabel('radians','Color',[1 0 1]);

subplot(5,1,5);
yinv=ifft(Yinv);
plot(yinv,'m')
axis ([1 512 -6 9]);
title(' y_i_n_v(t)');
    
figure(20+i)
HdB = -20*log2(abs(H))/N;
subplot(6,2,[1 2]);
plot(HdB);
axis([1 512 0 .3]);
title('Absolute value of the transfer function H(w)');
ylabel('-20*log(H)/N','Color',[1 0 1]);

subplot(6,1,2);
plot(Phase_H,'g');
axis([1 512 -5 5]);
title('Phase of Transfer Function H(w)');
ylabel('radians','Color',[1 0 1]);

Log_Ywin(i,:) = -20*log2(abs(Ywin))/N;
subplot(6,1,3);
plot(Log_Ywin(i,:),'r');
axis([1 512 -0.2 .3]);
title('Absolute Value of Transfer Function Y_w_i_n(w)');
ylabel('-20*log(H)/N','Color',[1 0 1]);

SNR = SNR/N;
subplot(6,1,4);

plot(phase(Ywin),'m');
axis([1 512 -5 5]);
title('Phase of Transfer Function Y_w_i_n(w)');
ylabel('radians','Color',[1 0 1]);

SNR = SNR/N;
subplot(6,1,5);
plot(abs(1./SNR),'m');
axis([1 512 -.05 .05])
axis tight;
sTitle=sprintf('Noise to signal ratio \\phi_n_/_o(w)');
hTitle=title(sTitle);
set(hTitle,'Color',[0 0 1]);
ylabel('1/N','Color',[1 0 1]);

subplot(6,1,6);
ywin=ifft(Ywin);
plot(ywin,'m')
axis ([1 512 -.2 .2]);
title(' y_w_i_n(t)');




% plot charachteristic of
figure(24+i)
subplot(5,2,[1 2])
plot(HdB);
axis([1 512 0 .3]);
title('Absolute value of the transfer function H(w)');
ylabel('-20*log(H)/N','Color',[1 0 1]);


Phase_H = phase(H_new);
subplot(5,1,2);
plot(Phase_H,'g');
axis([1 512 -5 5]);
title('Phase of Transfer Function H(w)');
ylabel('radians','Color',[1 0 1]);

Log_Yhom = -20*log2(abs(Yhom))/N;
subplot(5,1,3);
plot(Log_Yhom,'r');
axis([1 512 -0.2 .3]);
title('Absolute Value of Transfer Function Y_h_o_m(w)');
ylabel('-20*log(H)/N','Color',[1 0 1]);

subplot(5,1,4);
plot(phase(Yhom),'y');
axis([1 512 -5 5]);
title('Phase of Transfer Function Y_h_o_m(w)');
ylabel('radians','Color',[1 0 1]);

subplot(5,1,5);
yhom=ifft(Yhom);
plot(yhom,'b')
axis ([1 512 -.2 .2]);
title(' y_h_o_m(t)');
