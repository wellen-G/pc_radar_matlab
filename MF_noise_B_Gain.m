%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%���۷�����matlab����ȹ������ܶȵĸ�˹����������ӳ��ͬ�����µ��������ʱ仯��
%���ź�Ϊʱ��10us����30MHz/15MHz��LFM���壬������Ϊ80MHz����ӳ�źŴ������ѹ�����Ӱ�졣
% 2021/11/30
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
fs=8e7;                    %����Ƶ��Ϊ80MHz
Ts=1/fs;
B=15E6;                    %LFM�źŴ���Ϊ30/15MHz
T=10e-6;                   %LFM�ź�ʱ��Ϊ10us
N=1024;                    %FFT�ĵ���
SNR_Preset=-10;              %�����ź������
k=B/T;                     %��Ƶб��
A_xnoise=1;                  %������׼��
n=round(T*fs);             %���������
%
%    �������ѹ���棺BT      ��ѹ�������������źŵ�ʱ�������й�
%       Gain=10*log10(BT)
%       10*log10(300)=24.7712     10*log10(150)=21.7609
%    ���ı��źŵĴ�����ܸı���������ȣ��ı䲻����������
%    ��ѹ������������ʱ���й�
%    
%
%    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=linspace(-T/2,T/2,n);
s_Preset=exp(1j*pi*k*t.^2);          %����ʱ��Ϊ10us,����Ϊ30MHz�����Ե�Ƶ�ź�
h=conj(fliplr(s_Preset));            %ƥ���˲����ĳ弤��Ӧ
xnoise_Preset= A_xnoise*randn(size(s_Preset))+1i*A_xnoise*randn(size(s_Preset));  %�������źŵȳ��ȵ�����

if B<fs
    %%%%%     ��ư�����˹��ͨ�˲�    %%%%%%%%%%%%%%%%%%%%%%%%%
    Wc=B/fs;                                          %��ֹƵ��ΪB/2Hz
    [b,a]=butter(20,Wc);
    xnoise=filter(b,a,xnoise_Preset);    %���������е�ͨ�˲�
else
    xnoise=xnoise_Preset;
end


 
%%%%%%%%%%%%%%%%%%%%%%%%  �źŴ������������ȵ�Ӱ��  %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%    �����������ܶȲ���  %%%%%%%%%%%%%%%%%%%%%%   
[SNR_in,p_s,p_n]=SNR_caj(s_Preset,xnoise);
[s_mf_out,xnoise_mf_out]=FFT_MF(s_Preset,xnoise,h);
[SNR_out,p_s_max]=SNR_max_out(s_mf_out,xnoise_mf_out);
%%%%%%%%%%%%%%%%%%%%    �̶����������  %%%%%%%%%%%%%%
[s,s_add_noise,SNR_in,xnoise]=add_SNR_xnoise(s_Preset,SNR_Preset,xnoise);
h=conj(fliplr(s));  %ƥ���˲����ϵͳ����
[SNR_in_2,p_s,p_n]=SNR_caj(s,xnoise);
[s_mf_out,xnoise_mf_out]=FFT_MF(s,xnoise,h);
[SNR_out,p_s_max]=SNR_max_out(s_mf_out,xnoise_mf_out);

%%%%%%%%%%%%%%% ��ͼ %%%%%%%%%%%%%%%%%%%%%%%%%%%%5 
figure(1);   %  �������Ե�Ƶ�źŵ�ʱ��� Ƶ����
subplot(211);
plot(t,real(s_Preset)),title('LFM�ź�ʱ��'),xlabel('t/s'),ylabel('����');
S=fftshift(fft(s_Preset,N));
f=linspace(-fs/2,fs/2,N);
subplot(212);
plot(f,abs(S)),title('LFM�ź�Ƶ��'),xlabel('f/Hz'),ylabel('����');

figure(2);   %  �������Ե�Ƶ�źŵ�ʱ��� Ƶ����
subplot(211);
plot(t,real(xnoise)),title('����ʱ��'),xlabel('t/s'),ylabel('����');
NK=fftshift(fft(xnoise,N));
f=linspace(-fs/2,fs/2,length(NK));
subplot(212);
plot(f,abs(NK)),title('����Ƶ��'),xlabel('f/Hz'),ylabel('����');




function [s,s_add_noise,SNR_in,xnoise]=add_SNR_xnoise(s_Preset,SNR_Preset,xnoise)
    %  add_SNR_xnoise���ź��м���̶����������
    %  s_PresetΪ�������������ź�    SNR_PresetΪ�����  A_xnoiseΪ�����ı�׼��
    %  sΪ�ź�     s_add_noise ������������ź�   SNR_in���������    xnoiseΪ����
    Pn_add=sum(xnoise .*conj(xnoise ));                % ����������������� 
    Ps_add=sum(s_Preset.*conj(s_Preset));               % ���㵥λ�źŵ����� 
    A_1= sqrt((10^(SNR_Preset/10))*(Pn_add/Ps_add));     %�����Ӧ��������źŵķ���
    s=A_1.*s_Preset;                                    %�õ���Ӧ������µ��ź�
    s_add_noise=s+xnoise;                            % �����������ź�
    Ps_add_out=sum(s.*conj(s));
    SNR_in=10*log10(Ps_add_out/Pn_add);                   %���������
end


function   [s_mf_out,xnoise_mf_out]=FFT_MF(s_in,xnoise_in,h_in)
   % FFTʵ��ƥ���˲�
   
   xnoise_in=[xnoise_in,xnoise_in,xnoise_in,xnoise_in];   %���������ӳ�
  
   N1=length(s_in);
   N2=length(h_in);
   N_sum=N1+N2-1;
   v=ceil(log2(N_sum));
   N_s=2^v;
   SK=fft(s_in,N_s);
   HK=fft(h_in,N_s);
   SK_out=SK.*HK;          %LFM�ź���MF��FFTȻ�����
   s_mf_out=ifft(SK_out,N_s);
   
   N3=length(xnoise_in);
   N_n=2^ceil(log2(N1+N3-1));
   NK=fft(xnoise_in,N_n);
   HK_n=fft(h_in,N_n);
   NK_out=NK.*HK_n;
   xnoise_mf_out_temp=ifft(NK_out,N_n);
   xnoise_mf_out=xnoise_mf_out_temp(N_n/2-N_n/2*(N_s/N_n)+1:N_n/2+N_n/2*(N_s/N_n));
   
   
end

function  [SNR_out,p_s_max]=SNR_max_out(s,xnoise)
  %���������ȼ���ֵ����
  p_s_max=10*log10(max(s)*conj(max(s)));
  p_n=10*log10(sum(xnoise.*conj(xnoise))/length(xnoise));
  SNR_out =p_s_max- p_n;
end

function  [SNR_out,p_s,p_n]=SNR_caj(s,xnoise)
  %������ȼ�ƽ������
  p_s=10*log10(sum(s.*conj(s))/length(s));
  p_n=10*log10(sum(xnoise.*conj(xnoise))/length(xnoise));
  SNR_out =p_s- p_n;
end




