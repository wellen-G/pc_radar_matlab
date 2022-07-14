%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ������������5/10MHz������Ƶ���ź�Ƶ����500KHz����������5MHz���������FFT�����档
%   SNR�ɵ�   FFT�����ɵ�
% 2021/11/30
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
clc
N=1024;        % FFT�ĵ���Ϊ   N=1024
fs=10E6;        % ������       fs=5MHz   
Ts=1/fs;
%     ��ʱ�ֱ���Ϊfs/N= 4.8828e+03    4.8KHz
B_noise=10e6;   %��������
fd=1.25e6;     % �ź�Ƶ��     fd=1.25MHz
% fd=1.25e6+D(i_avr)*fs/N;     % �ź�Ƶ��     fd=1.25MHz
% C=0.5;
% fd=1.25e6+C*fs/N;
SNR_Preset=-10;         %��ʼ����ȣ���dBΪ��λ
A_xnoise=1;                  %������׼��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=0:Ts:N*Ts-Ts;
n=length(t); %������ĸ���
s_Preset=exp(1i*(2*pi*fd*t+pi/6));      %�ȶ����ֵΪ1��Ƶ���ź�
% s_Preset=exp(1i*(2*pi*fd*t+pi/6))+exp(1i*(2*pi*1.25e6*t+pi/6)); 
xnoise_Preset= A_xnoise*randn(size(s_Preset))+1i*A_xnoise*randn(size(s_Preset));  %�������źŵȳ��ȵ�����
if B_noise<fs
    %%%%%     ��ư�����˹��ͨ�˲�    %%%%%%%%%%%%%%%%%%%%%%%%%
    Wc=B_noise/fs;                                          %��ֹƵ��ΪB/2Hz
    [b,a]=butter(20,Wc);
    xnoise=filter(b,a,xnoise_Preset);    %���������е�ͨ�˲�
else
    xnoise=xnoise_Preset;
end

%%%%%%%%%%%%%%  �źŲ���ʱ��������ı仯����������ȵ�Ӱ��  %%%%%%%%%%%%%%55
[SNR_in,p_s,p_n]=SNR_caj_time(s_Preset,xnoise);   
% s=s_Preset;
[s,s_add_noise,SNR_in,xnoise]=add_SNR_xnoise(s_Preset,SNR_Preset,xnoise);
%%%%%%%%%%%%%%%%%%%%  ��֤�����غ�  %%%%%%%%%%%%%%%%%%%%%%%
[SNR_in,p_s,p_n]=SNR_caj_time(s,xnoise);   
SK=fftshift(fft(s,N));
NK=fftshift(fft(xnoise,N)); 
% XK=fftshift(fft(s_add_noise,N)); 
[SNR_K_out,P_SK,P_NK]=SNR_caj_fre(SK,NK); 

[SNR_FFT_out,p_s_max,p_NK_2]=SNR_max_out(SK,NK,fs,B_noise);  % FFT��ķ�ֵ�����
Gain_FFT=SNR_FFT_out-SNR_in;
%%%%%%%%%%%%%%%%%% Ƶ��λ�ú�FFT���������ȵĹ�ϵ %%%%%%%%%%%%%%%%%%%%%%%%%%
%  ����avr_SNR������ִ��

%%%%%%%%%%%%%%%%%%%%% ��ͬ���������������ȵ�Ӱ��  %%%%%%%%%%%%%%%%%   


 w(:,1)=hann(length(t));      %����������
 w(:,2)=chebwin(length(t));   %�б�ѩ�򴰺���
%   wvtool(w(:,2))
 w(:,3)=triang(length(t));    %���Ǵ�����
 w(:,4)=hamming(length(t));   %����������
 w(:,5)=blackman(length(t));  %��������������
 w(:,6)=boxcar(length(t));     %���δ��������൱��δ�Ӵ�)
  
w_N=length(w(1,:));                    %���봰����������
xnoise_addw=zeros(w_N,length(s));      %�Ӵ��������
s_addw=zeros(w_N,length(s));           %�Ӵ�����ź�
SK_addw=zeros(w_N,N);                  %FFT����ź�
NK_addw=zeros(w_N,N);                  %FFT�������
SNR_in_addw=zeros(1,w_N);              %FFTǰ�ķ�ֵ�����
SNR_out_addw=zeros(1,w_N);             %FFT��ķ�ֵ�����
Gain_in_addw=zeros(1,w_N);             %FFT��ķ�ֵ���������
P_s_max_addw=zeros(1,w_N);             %FFT����źŷ�ֵ����

for i=1:w_N
  xnoise_addw(i,:)=xnoise.*(w(:,i)');    %�����Ӵ�
  s_addw(i,:)=s.*w(:,i)';                %�źżӴ�
  [SNR_in_addw(1,i),~,~]=SNR_caj_time(s_addw(i,:),xnoise_addw(i,:));   
  SK_addw(i,:)=fftshift(fft(s_addw(i,:),N));
  NK_addw(i,:)=fftshift(fft(xnoise_addw(i,:),N)); 
  [SNR_out_addw(1,i),~,~]=SNR_max_out( SK_addw(i,:),NK_addw(i,:),fs,B_noise);  % FFT��ķ�ֵ�����
  Gain_in_addw(1,i)=SNR_out_addw(1,i)-SNR_in_addw(1,i);          %FFT��ķ�ֵ���������
end


%%%%%%%%%%%%%%  ��ͼ %%%%%%%%%%%%%%%%%%%%%%%
% f=linspace(-fs/2,fs/2-fs/length(SK),length(SK));
% figure(1)   %%%%%%%%%  �ź�  %%%%%%%%%%
% % subplot(311)
% % plot(t,real(s)), xlabel('t/s'),ylabel('����'),title('�ź�ʱ����'),grid on;
% 
% subplot(211)
% plot(t,real(s),'linewidth',1), xlabel('t/s'),ylabel('����'), xlim([0 ,200*Ts]),title('�ź�ʱ����'),grid on;
% 
% subplot(212)
% plot(f,abs(SK),'linewidth',1), xlabel('f/Hz'),ylabel('����'),title('�ź�Ƶ������'),grid on;
% 
% % figure(2)        %%%%%%%%%%%  ����  %%%%%%%%%%
% % subplot(211)
% % plot(t,abs(xnoise)), xlabel('t/s'),ylabel('����'),title('����ʱ����'),grid on;
% % 
% % subplot(212)
% % 
% % plot(f,abs(NK),'linewidth',1), xlabel('f/Hz'),ylabel('����'),title('����Ƶ������'),grid on;
% 
% figure(3)        %%%%%%%%%%%  �ź�+����  %%%%%%%%%%
% % subplot(311)
% % plot(t,real(s_add_noise)), xlabel('t/s'),ylabel('����'),title('�ź�+����ʱ����'),grid on;
% subplot(211)
% plot(t,real(s_add_noise),'linewidth',1), xlabel('t/s'),ylabel('����'),xlim([0 ,200*Ts]),title('�ź�+����ʱ����'),grid on;
% subplot(212)
% plot(f,abs(XK),'linewidth',1), xlabel('f/Hz'),ylabel('����'),title('�ź�+����Ƶ������'),grid on;
% 
figure(4)  
subplot(211)
plot(f,abs(SK),'linewidth',1.2), xlabel('f/Hz'),ylabel('����'),title('�ź�Ƶ������'),grid on;
subplot(212)
stem(f,abs(SK) ,'linewidth',1.2), xlabel('f/Hz'),ylabel('����'),xlim([1.2e6 ,1.3e6]),title('�ź�Ƶ������(С������)'),grid on;



function [s,s_add_noise,SNR_in,xnoise]=add_SNR_xnoise(s_Preset,SNR_Preset,xnoise)
    %  add_SNR_xnoise���ź��м���̶����������
    %  s_PresetΪ�������������ź�    SNR_PresetΪ�����  A_xnoiseΪ�����ı�׼��
    %  sΪ�ź�     s_add_noise ������������ź�   SNR_in���������    xnoiseΪ���� 
    Pn_add=sum(xnoise .*conj(xnoise ));                  % ������������Ĺ��� 
    Ps_add=sum(s_Preset.*conj(s_Preset));                % ���㵥λ�źŵĹ��� 
    A_1= sqrt((10^(SNR_Preset/10))*(Pn_add/Ps_add));     %�����Ӧ��������źŵķ���
    s=A_1.*s_Preset;                                     %�õ���Ӧ������µ��ź�
    s_add_noise=s+xnoise;                                % �����������ź�
    Ps_add_out=sum(s.*conj(s));
    SNR_in=10*log10(Ps_add_out/Pn_add);                   %���������
end
function  [SNR_out,p_s,p_n]=SNR_caj_fre(s,xnoise)
  %������ȼ�ƽ������
  p_s=10*log10(sum(s.*conj(s))/length(s));
  p_n=10*log10(sum(xnoise.*conj(xnoise))/length(xnoise));
  SNR_out =p_s- p_n;
end
function  [SNR_out,p_s,p_n]=SNR_caj_time(s,xnoise)
  %������ȼ�����
  p_s=10*log10(sum(s.*conj(s)));
  p_n=10*log10(sum(xnoise.*conj(xnoise)));
  SNR_out =p_s- p_n;
end
function  [SNR_out,p_s_max,p_n]=SNR_max_out(s,xnoise,fs,B_n)
  %���������ȼ���ֵ����
  %fsΪ������  B_nΪ��������
  n=length(xnoise);
  n1=round(n/2-n/2*(B_n/fs))+1;
  n2=round(n/2+n/2*(B_n/fs));
  p_s_max=10*log10(max(s)*conj(max(s)));
  p_n=10*log10(sum(xnoise(n1:n2).*conj(xnoise(n1:n2)))/(n2-n1+1));
  SNR_out =p_s_max- p_n;
end