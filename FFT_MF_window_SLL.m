%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ƥ���˲�����FFT-IFFT��ʵ��
%   ���ź�Ϊʱ��10us ����30MHz��LFM���壬������Ϊ80MHz����ѧ�����ͷ��档
%   �����Ӵ�����԰�ı仯 ����ȵ���ʧ
%   2021/12/1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
close all
fs=8e7;               %����Ƶ��Ϊ80MHz        ����Ƶ��Ϊ10MHzʱ��������ӽ�����ֵ
Ts=1/fs;
B=30E6;              %LFM�źŴ���Ϊ30MHz
T=10e-6;             %LFM�ź�ʱ��Ϊ10us
% N=2048;            %FFT�ĵ���
SNR=-10;
k=B/T;               %��Ƶб��
A_xnoise=1;
n=round(T*fs);       %���������
N=1024;                    %FFT�ĵ���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
t=linspace(-T/2,T/2,n);
s_Preset=exp(1j*pi*k*t.^2);       %����ʱ��Ϊ10us,����Ϊ30MHz�����Ե�Ƶ�ź�
xnoise_Preset= A_xnoise*randn(size(s_Preset))+1i*A_xnoise*randn(size(s_Preset));  %�������źŵȳ��ȵ�����

if B<fs
    %%%%%     ��ư�����˹��ͨ�˲�    %%%%%%%%%%%%%%%%%%%%%%%%%
    Wc=B/fs;                                          %��ֹƵ��ΪB/2Hz
    [b,a]=butter(20,Wc);
    xnoise=filter(b,a,xnoise_Preset);    %���������е�ͨ�˲�
else
    xnoise=xnoise_Preset;
end
 %xnoise=xnoise_Preset;
%plot(t,s_Preset)

% s=s_Preset;SNR_in=10*log10(sum(s.*conj(s))/sum(xnoise.*conj(xnoise)));
% p_ss=10*log10(sum(s.*conj(s))/length(s));p_nn=10*log10(sum(xnoise.*conj(xnoise))/length(xnoise));
[s,s_add_noise,SNR_in,xnoise]=add_SNR_xnoise(s_Preset,SNR,xnoise);   %�õ�����������µ��ź��Լ�����
h=conj(fliplr(s));      %ƥ���˲���ĳ弤��Ӧ

[s_out,xnoise_out]=FFT_MF(s,xnoise,h);               %FFTʵ��ƥ���˲� 
[SNR_out,P_s_max,p_n_out]=SNR_max_out(s_out,xnoise_out);     %������ѹ��ķ�ֵ����Ⱥ��źŷ�ֵ����
Gain_MF=SNR_out-SNR_in;
%%%%%%%%%%%%%%%%%%%%% ���ǻز��뷢���źŲ��غ� %%%%%%%%%%%%%%%%%
N_return=40;
for i=1:N_return
    
    t_1=linspace(-T/2-Ts*(i/(N_return)),T/2-Ts*(i/(N_return)),n);
    s_reture_2(i,:)=exp(1j*pi*k*t_1.^2); 
    PS_reture_in_2(i)=10*log10(sum(s_reture_2(i,:).*conj(s_reture_2(i,:)))/length(s_reture_2(i,:)));
    SNR_in_reture_2(i)=PS_reture_in_2(i)-10*log10(sum(xnoise.*conj(xnoise))/length(xnoise));
%     [s_reture_2(i,:),s_add_noise_reture_2(i,:),SNR_in_reture_2(i),xnoise]=add_SNR_xnoise(s_return_2(i,:),SNR,xnoise);   %�õ�����������µ��ź��Լ�����
    [s_out_reture_2(i,:),xnoise_out]=FFT_MF(s_reture_2(i,:),xnoise,h);               %FFTʵ��ƥ���˲� 
    [SNR_out_reture_2(i),P_s_max_reture_2(i)]=SNR_max_out(s_out_reture_2(i,:),xnoise_out);     %������ѹ��ķ�ֵ����Ⱥ��źŷ�ֵ����
    Gain_MF_reture(i)=SNR_out_reture_2(i)-SNR_in_reture_2(i);
    
end

figure(6);   % �ز��뷢���źŲ��غ�ʱ
t_n=linspace(0,1,N_return);
plot(t_n,Gain_MF_reture,'linewidth',1),title('�ز��뷢���źŲ��غ�ʱ,��ѹ�����������Ĺ�ϵ'),xlabel('��������ĸ���'),ylabel('����/dB');
clear i N_return;
    
t_1=linspace(-T/2-Ts/2,T/2-Ts/2,n);
s_return=exp(1j*pi*k*t_1.^2); 
[s_reture,s_add_noise_reture,SNR_in_reture,xnoise]=add_SNR_xnoise(s_return,SNR,xnoise);   %�õ�����������µ��ź��Լ�����
[s_out_reture,xnoise_out]=FFT_MF(s_reture,xnoise,h);               %FFTʵ��ƥ���˲� 
[SNR_out_reture,P_s_max_reture]=SNR_max_out(s_out_reture,xnoise_out);     %������ѹ��ķ�ֵ����Ⱥ��źŷ�ֵ����
Gain_MF_reture=SNR_out_reture-SNR_in_reture;


%%%%%%%%%%%%%%%%%%%%%% ���мӴ�  %%%%%%%%%%%%%%%%%%%%%%%%%

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
s_out_addw=zeros(w_N,2^round(log2(2*length(s)-1)));       %��ѹ����ź�
xnoise_out_addw=zeros(w_N,2^round(log2(2*length(s)-1)));   %��ѹ�������
SNR_out_addw=zeros(1,w_N);             %��ѹ��ķ�ֵ�����
P_s_max_addw=zeros(1,w_N);             %��ѹ����źŷ�ֵ����
SLL_addw=zeros(1,w_N);                 %��ѹ����ź��԰��ƽ
SLL_addw_interp=zeros(1,w_N);           %��ѹ����ź��԰��ƽ_����ֵ��
for i=1:w_N
  xnoise_addw(i,:)=xnoise.*(w(:,i)');
  s_addw(i,:)=s.*w(:,i)';           %���뺺����
  h_addw=h.*w(:,i)';
  [s_out_addw(i,:),xnoise_out_addw(i,:)]=FFT_MF(s_addw(i,:),xnoise_addw(i,:),h);      %FFTʵ��ƥ���˲� 
  [SNR_out_addw(i),P_s_max_addw(i)]=SNR_max_out(s_out_addw(i,:),xnoise_out_addw(i,:));     % �õ����������Լ��źŷ�ֵ����
  [SLL_addw(i)]=SLL_CAL(s_out_addw(i,:));     %�����԰��ƽ
  
end


%%%%%%%%%%%%%%%%%%   ����ֵ���� %%%%%%%%%%%%%%%%%%%%
N_interp=10;      %�ڲ������   
s_out_addw_interp=zeros(w_N,2^round(log2(2*length(s)-1))*N_interp);       %��ѹ����ź�
for i=1:w_N
    s_out_addw_interp(i,:)=interp(s_out_addw(i,:),N_interp);
    [SLL_addw_interp(i)]=SLL_CAL(s_out_addw_interp(i,:));        %��ֵ����԰��ƽ
end
 
Ts_interp=Ts/N_interp;

[s_out_addw_peak_1,~,tmax_1,~]=Peak(s_out_addw_interp(1,:),Ts_interp);     %��ѹ����źż�ֵ��'
[s_out_addw_peak_2,~,tmax_2,~]=Peak(s_out_addw_interp(2,:),Ts_interp);     %��ѹ����źż�ֵ��'
[s_out_addw_peak_3,~,tmax_3,~]=Peak(s_out_addw_interp(3,:),Ts_interp);     %��ѹ����źż�ֵ��'
[s_out_addw_peak_4,~,tmax_4,~]=Peak(s_out_addw_interp(4,:),Ts_interp);     %��ѹ����źż�ֵ��'
[s_out_addw_peak_5,~,tmax_5,~]=Peak(s_out_addw_interp(5,:),Ts_interp);     %��ѹ����źż�ֵ��'
[s_out_addw_peak_6,~,tmax_6,~]=Peak(s_out_addw_interp(6,:),Ts_interp);     %��ѹ����źż�ֵ��'

%  
% figure(1)
% plot(t,real(h.*w(:,2)'))

%%%%%%%%%%%%%%%%%%%%%5 ��ͼ  %%%%%%%%%%%%%%%%%%%%%%%%%%%%55
figure(1);   %  �������Ե�Ƶ�źŵ�ʱ��� Ƶ����
subplot(211);
plot(t,real(s_Preset)),title('LFM�ź�ʱ��'),xlabel('t/s'),ylabel('����');
S=fftshift(fft(s_Preset,N));
f=linspace(-fs/2,fs/2,N);
subplot(212);
plot(f,abs(S)),title('LFM�ź�Ƶ��'),xlabel('f/Hz'),ylabel('����');

figure(4);   %  �������Ե�Ƶ�źŵ�ʱ��� Ƶ����
subplot(211);
plot(t,real(xnoise)),title('����ʱ��'),xlabel('t/s'),ylabel('����');
NK=fftshift(fft(xnoise,N));
f=linspace(-fs/2,fs/2,length(NK));
subplot(212);
plot(f,abs(NK)),title('����Ƶ��'),xlabel('f/Hz'),ylabel('����');

figure(5);   % �ز��뷢���źŲ��غ�ʱ
hold on
plot(t,real(s_Preset),'-*','linewidth',1),title('LFM�ź�ʱ��'),xlabel('t/s'),ylabel('����');
plot(t,real(s_return),'-+','linewidth',1),title('LFM�ź�ʱ��'),xlabel('t/s'),ylabel('����'),ylim([-1.5,1.5]);
legend('�����ź�','�ز��ź�')
hold off


figure(3)
hold on
t_2=Ts:Ts:(length(s_out))*Ts;
plot(t_2,abs(s_out)),xlim([0.9E-5,1.1E-5]),title('ƥ���˲����ź�ʱ���Σ�FFT��')
plot(t_2,abs(xnoise_out)),xlim([0.9E-5,1.1E-5])
plot(t_2,abs(xnoise_out+s_out)),xlim([0.9E-5,1.1E-5])
hold off

figure(2)
t_interp=Ts_interp:Ts_interp:(length(s_out_addw_interp(1,:)))*Ts_interp;
m_1=0.95E-5;m_2=1.05E-5;
subplot(231)
hold on
plot(t_interp,20*log10(abs(s_out_addw_interp(1,:)))),xlim([m_1,m_2]),title('ƥ���˲����ź�ʱ���Σ���������')
plot(tmax_1,20*log10(abs(s_out_addw_peak_1)),'r+') 
hold off
subplot(232)
hold on
plot(t_interp,20*log10(abs(s_out_addw_interp(2,:)))),xlim([m_1,m_2]),title('ƥ���˲����ź�ʱ���Σ��б�ѩ�򴰣�')
plot(tmax_2,20*log10(abs(s_out_addw_peak_2)),'r+')
hold off
subplot(233)
hold on
plot(t_interp,20*log10(abs(s_out_addw_interp(3,:)))),xlim([m_1,m_2]),title('ƥ���˲����ź�ʱ���Σ����Ǵ���')
plot(tmax_3,20*log10(abs(s_out_addw_peak_3)),'r+')
hold off
subplot(234)
hold on
plot(t_interp,20*log10(abs(s_out_addw_interp(4,:)))),xlim([m_1,m_2]),title('ƥ���˲����ź�ʱ���Σ���������')
plot(tmax_4,20*log10(abs(s_out_addw_peak_4)),'r+')
hold off
subplot(235)
hold on
plot(t_interp,20*log10(abs(s_out_addw_interp(5,:)))),xlim([m_1,m_2]),title('ƥ���˲����ź�ʱ���Σ�������������')
plot(tmax_5,20*log10(abs(s_out_addw_peak_5)),'r+')
hold off
subplot(236)
hold on
plot(t_interp,20*log10(abs(s_out_addw_interp(6,:)))),xlim([m_1,m_2]),title('ƥ���˲����ź�ʱ���Σ�δ�Ӵ���')
plot(tmax_6,20*log10(abs(s_out_addw_peak_6)),'r+')
hold off

 [SLL_addw_interp(5)]=SLL_CAL(s_out_addw_interp(5,:));  
 
function [s,s_add_noise,SNR_in,xnoise]=add_SNR_xnoise(s_Preset,SNR_Preset,xnoise)
    %  add_SNR_xnoise���ź��м���̶����������
    %  s_PresetΪ�������������ź�    SNR_PresetΪ�����  A_xnoiseΪ�����ı�׼��
    %  sΪ�ź�     s_add_noise ������������ź�   SNR_in���������    xnoiseΪ���� 
    Pn_add=sum(xnoise .*conj(xnoise));                % ����������������� 
    Ps_add=sum(s_Preset.*conj(s_Preset));             % ���㵥λ�źŵ����� 
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
   N3=length(xnoise_in);
   N_s=2^ceil(log2(N1+N2-1));
   SK=fft(s_in,N_s);
   HK=fft(h_in,N_s);
   SK_out=SK.*HK;          %LFM�ź���MF��FFTȻ�����
   s_mf_out=ifft(SK_out,N_s);
  
   N_n=2^ceil(log2(N3+N2-1));
   NK=fft(xnoise_in,N_n);
   HK_n=fft(h_in,N_n);
   NK_out=NK.*HK_n; 
   xnoise_mf_out_temp=ifft(NK_out,N_n);
   xnoise_mf_out=xnoise_mf_out_temp(N_n/2-N_n/2*(N_s/N_n)+1:N_n/2+N_n/2*(N_s/N_n));
   
end

function  [SNR_out,p_s_max,p_n]=SNR_max_out(s,xnoise)
  %���������ȼ���ֵ����
  p_s_max=10*log10(max(s)*conj(max(s)));
  p_n=10*log10(sum(xnoise.*conj(xnoise))/length(xnoise));
  SNR_out =p_s_max- p_n;
end


function  [vmax,vmin,tmax,tmin]=Peak(signal,Ts)
   %  �����еļ�ֵ
    t=Ts:Ts:(length(signal))*Ts;
    Lmax=diff(sign(diff(abs(signal))))==-2;
    Lmin=diff(sign(diff(abs(signal))))==2;
    Lmax=[false,Lmax,false];
    Lmin=[false,Lmin,false];
    tmax=t(Lmax);
    tmin=t(Lmin);
    vmax=signal(Lmax);
    vmin=signal(Lmin);
end

function [SLL]=SLL_CAL(signal)
    %���źŵ��԰��ƽ
   
%     signal=signal(round(11/30*length(signal)):round(19/30*length(signal)));
    Lmax=diff(sign(diff(abs(signal))))==-2;    %�õ���ֵ��
    Lmax=[false,Lmax,false];
    vmax=signal(Lmax);
    [S_max,S_max_index]=max(vmax);             %��������ƽ��λ��
    S_sec_sum=[vmax(S_max_index-5:S_max_index-1),vmax(S_max_index+1:S_max_index+5)];  %�õ�������Χ���԰�
    S_sec=max(S_sec_sum);               %�õ�����԰�
%     Signal_sort=sort(abs(vmax));    %����������
%     S_max= Signal_sort(length(Signal_sort));
%     S_sec=Signal_sort(length(Signal_sort)-1);
    SLL=20*log10(abs(S_sec)/abs(S_max));
 end
