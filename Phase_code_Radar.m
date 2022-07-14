clc 
clear 
close all 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m_length=127;                     %�볤Ϊ127
fm=72e6;                          %��ƵΪ72MHz
fc=10e9;                          %�ź���ƵΪ 10GHz
fs=72e6*4;                         %������Ϊ200MHz
Ts=1/fs;
Duty_Ratio= 0.10;                  %ռ�ձ�Ϊ10%
PRT=m_length*(1/fm)/Duty_Ratio;    %�����ظ�����
PRF=1/PRT;                        %�����ظ�Ƶ��
n_add=round(100e-3/(m_length/fm)*Duty_Ratio);        %����ۼƸ���
n_add=2^round(log2(n_add));

 n_add=256;
CPI=n_add*PRT;                       %��δ���ʱ��
SNR_Preset=0;                       %���������
A_xnoise =1;                            %������׼��
N_PRT=m_length*(fs/fm)/Duty_Ratio;    %һ�������ڵĲ�������
T_width=m_length/fm;              %�źŵ�ʱ��
B=fm;
%%%%%%%%%%%%%%%%%%%%%%  ����m����   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ��˹����ź� 
u = idinput(m_length, 'prbs',[0,1],[-1,1])';   %����m����
%  figure(12) ,stairs(u,'linewidth',1) ,title('m����')
[c,lags] = xcorr(u);    % ��һ����m���е�����غ���
%  figure(11),  plot(lags,c,'linewidth',1),title('m���������')   %����m���е������
 m_sum=[];
 z=zeros(m_length*10*(1-Duty_Ratio),1);
 for i=1:n_add  
    m_sum=[m_sum,idinput(m_length, 'prbs',[0,1],[-1,1])'];%,z'];   
 end
[c_sum,lags_sum] = xcorr(u ,m_sum);    % ��һ����m���е�ѭ������غ���
 % figure(10) , plot(c_sum,'linewidth',1),title('m���е�ѭ�������')     %����m���е�ѭ�������
%%%%%%%%%%%%%%%%%%%%%%%% ����α�����λ��������ź� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

Phase_code_signal=Phase_code_modulation( u,fm,fc,fs,Duty_Ratio,n_add) ;      %���������ź�
%    plot( real(Phase_code_signal)),xlim([0,2e5]),ylim([-2,2]),ylim([-2,2]), title('α�����λ�����ź�') 
s_length=length(Phase_code_signal);                        %�����źŵĳ���
xnoise_Preset= A_xnoise*randn(size(Phase_code_signal))+1i*A_xnoise*randn(size(Phase_code_signal));  
 [signal_demodulation]=demodulation(Phase_code_signal,fc,fs);     %�Է����źŽ�� 
  h=conj(fliplr(signal_demodulation(1:m_length*round(fs/fm))));   %����ƥ���˲������Գ弤��Ӧ 
 

 SK_mf=fftshift(fft(signal_demodulation));
 f_s=linspace(-fs/2,fs/2,length(Phase_code_signal));
 % figure(1),plot(f_s,abs(SK_mf));

 %%%%%%%%% ��֤��ѹ������  %%%%%%%%%%%%%%%%%%%%%%%
   %%%%%     ��ư�����˹��ͨ�˲�    %%%%%%
  Wc=B/fs;                                          %��ֹƵ��ΪB/2Hz
  N_filter=20;     % �˲����Ľ���
  [b_filter,a_filter]=butter(N_filter,Wc);

  xnoise=filter(b_filter,a_filter,xnoise_Preset(1:N_PRT/10));    %���������е�ͨ�˲�
%   s=filter(b_filter,a_filter,signal_demodulation(1:N_PRT/10));     %���źŽ��е�ͨ�˲�
  [ signal_filter_return]=filter_compose(signal_demodulation,a_filter,b_filter,N_PRT,Duty_Ratio,n_add);
  [ xnoise_filter_return]=filter_compose(xnoise_Preset,a_filter,b_filter,N_PRT,Duty_Ratio,n_add);
  
  %%%%%%%%%%%%%%  ������ѹ������  %%%%%%%%%%%%%%%
  %  �Ե�һ������Ϊ��  ������ʱ   �޶�����ЧӦ
  [s_snr,s_add_noise,SNR_in,xnoise_snr]=add_SNR( signal_filter_return(1:N_PRT/10),SNR_Preset,xnoise_filter_return(1:N_PRT/10));
  
  [signal_mf_out,noise_mf_out]=FFT_MF(s_snr(1:N_PRT/10),xnoise_snr(1:N_PRT/10),h);
  [SNR_in_MF,~,~]=SNR_caj_time(s_snr,  xnoise_snr);
  [SNR_out_MF,~,~]=SNR_max_out(signal_mf_out,noise_mf_out);
  Gain_MF=SNR_out_MF-SNR_in_MF;
 
  NK_mf=fftshift(fft(xnoise));
  f_n=linspace(-fs/2,fs/2,length(xnoise));
 % figure(2),plot(f_n,abs(NK_mf));

%%%%%%%%%%%%%%%%%%%%%%%%%%%  �����ز��ź� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�ٶȷֱ���  
v_distinguish=3e8*PRF/2/n_add/fc;
%����ֱ���  3e8/fm/2
R_distinguish=3e8/fm/2;
R_max=(3e8/2)*PRT;     %���ģ������
v_max=3e8*PRF/2/fc;    %���ģ���ٶ�

R_target=[1000,1000];               %Ŀ��ľ���    R_max=(3e8/2)*PRT*(1-Duty_Ratio)=2381 R_min=(3e8/2)*PRT*Duty_Ratio=264
v_target=[400,435];                          %Ŀ����ٶ�    v_max=3e8*PRF/2/fc=850
A_target=[1,0.04];                        %Ŀ��ز��źŵķ�ֵ˥��

%Ŀ��һ���״�ز�
Phase_code_signal_return_1=Phase_code_signal_return_modulation(Phase_code_signal,fs,fc,Duty_Ratio,R_target(1),v_target(1),A_target(1),n_add);
%Ŀ������״�ز�
Phase_code_signal_return_2=Phase_code_signal_return_modulation(Phase_code_signal,fs,fc,Duty_Ratio,R_target(2),v_target(2),A_target(2),n_add);

%�Իز��źŽ��
 [Phase_code_signal_return_1]=demodulation(Phase_code_signal_return_1,fc,fs);
 [Phase_code_signal_return_2]=demodulation(Phase_code_signal_return_2,fc,fs);
 xnoise=filter(b_filter,a_filter,xnoise_Preset);    %���������е�ͨ�˲�
 % plot( real(Phase_code_signal_return_1)),xlim([0,2e5]),ylim([-2,2]), title('�޶�����ЧӦʱ�Ļز��ź�')
 
 % plot(abs(fft(Phase_code_signal_return_2)))
%����̶�����Ⱥ�Ļز��ź�
%�ϳɻز�+����
[Phase_code_signal_return_sum,Phase_code_signal_return_sum_add_noise,SNR_in_sum_return,xnoise_sum_return]=add_SNR_xnoise_twotarget(Phase_code_signal_return_1,Phase_code_signal_return_2,SNR_Preset,A_xnoise,m_length*round(fs/fm));
%Ŀ��һ�ز�+����
[Phase_code_signal_return_1 ,Phase_code_signal_return_1_noise, SNR_in_return(1),xnoise_return_1]=add_SNR_xnoise(Phase_code_signal_return_1,SNR_Preset,xnoise,m_length*round(fs/fm),N_filter);
%Ŀ����ز�+����
% [Phase_code_signal_return_2 ,Phase_code_signal_return_2_noise, SNR_in_return(2),xnoise_return_2]=add_SNR_xnoise(Phase_code_signal_return_2,SNR_Preset,A_xnoise,m_length*round(fs/fm));
 
% plot( real(Phase_code_signal_return_1_noise)), title('����������Ļز��ź�'),xlim([0,2e5]),ylim([-2,2])
% plot( real(Phase_code_signal_return_1)),xlim([0,2e5]),ylim([-2,2]), title('�޶�����ЧӦʱ�Ļز��ź�')

 SK_in=fftshift(fft(Phase_code_signal_return_1));
 SK_MF_in=fftshift(fft(xnoise_return_1));
 %��ѹǰ��Ĵ���
%  figure(8), subplot(211), plot(linspace(-fs/2,fs/2,length(SK_in)),abs(SK_in)),title('�ź�Ƶ��'),subplot(212), plot(linspace(-fs/2,fs/2,length( SK_MF_in)), abs(SK_MF_in)),title('����Ƶ��');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   ������ѹ ������ѹ  ʱ��  ����  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 [Phase_code_signal_return_1_mf_out,Phase_code_noise_return_1_mf_out]=FFT_MF(Phase_code_signal_return_1, xnoise_return_1,h);
%  [Phase_code_signal_return_2_mf_out,Phase_code_noise_return_2_mf_out]=FFT_MF(Phase_code_signal_return_2,xnoise_return_2,h);
 [Phase_code_signal_return_sum_mf_out,Phase_code_noise_return_sum_mf_out]=FFT_MF(Phase_code_signal_return_sum,xnoise_sum_return,h);
 
  % plot(real(Phase_code_signal_return_1_mf_out)), title('��ѹ����ź�'),xlim([0,2e5]),ylim([-2,2])
 
 %////////////////////  ��ѹ���ʱ��    ȡǰ�������ʱ�����  /////////////////////////////
 signal_mf_width_analysis=Phase_code_signal_return_1_mf_out(1:N_PRT*5);
 place_max=zeros(1,5);  %���ֵ����λ��
 place_4dB=zeros(1,5);  %�½�4dB����λ��

 N_interp=30;
%  s_interp = interp1( 1:1:N_PRT*5,signal_mf_width_analysis, 1/N_interp:1/N_interp:N_PRT*5,'spline')/max(signal_mf_width_analysis);
s_interp=interp(signal_mf_width_analysis, N_interp)/max(signal_mf_width_analysis);  %���в�ֵ����һ��
  for i=1:5
     [~,place_max(i)]=max(s_interp((i-1)*N_PRT* N_interp+1:i*N_PRT* N_interp));
     place_max(i)=place_max(i)+(i-1)*N_PRT *N_interp;
 end
  for i=1:5   %�ҳ��½�4dB����λ��
     [~,place_4dB(i)]=min(abs(10*log10(abs(s_interp((i-1)*N_PRT* N_interp+1:i*N_PRT* N_interp)))+4));
     place_4dB(i)=place_4dB(i)+(i-1)*N_PRT* N_interp;
  end
 T_mf=sum(abs(place_max- place_4dB))*2/fs/N_interp/5;   %��ѹ���ʱ��
 D_MF=T_width/T_mf;                                    %ѹ����
 
 t_interp=0:1/fs/N_interp:(length(s_interp)-1)/fs/N_interp;
%��ѹ���ʱ��
%  figure(9), plot(t_interp, 10*log10(abs(s_interp))),  xlim([(place_max(1)-150)*1/fs/N_interp,(place_max(1)+150)*1/fs/N_interp])
%//////////////////////// ��ѹ��Ĵ���/////////////////////////

 SK=fftshift(fft(Phase_code_signal_return_1));
 SK_MF=fftshift(fft( Phase_code_signal_return_1_mf_out));
 %��ѹǰ��Ĵ���
%  figure(8), subplot(211), plot(linspace(-fs/2,fs/2,length(SK)),abs(SK)),title('��ѹǰƵ��'),subplot(212), plot(linspace(-fs/2,fs/2,length( SK_MF)), abs(SK_MF)),title('��ѹ��Ƶ��');
%  ������ѹ�������Ϊ��Ƶfm
 %%%%%%%%%%%%%%%%%%%%%%%%%%%  ����������  %%%%%%%%%%%%%%%%%%%%%%%%%%5 
 n_y=n_add;                                %����ۼƸ���
 n_x=m_length*(fs/fm)/Duty_Ratio;          %�������ڲ�������
 signal_return_1_array=zeros(n_x,n_y);
 noise_return_1_array=zeros(n_x,n_y);
%  signal_return_2_array=zeros(n_y,n_x);
 signal_return_sum_array=zeros(n_x,n_y);

for i=1:n_y
    for j=1:n_x
        signal_return_1_array(j,i)=Phase_code_signal_return_1_mf_out((i-1)*n_x+j);
        noise_return_1_array(j,i)=Phase_code_noise_return_1_mf_out((i-1)*n_x+j);
        signal_return_sum_array(j,i)=Phase_code_signal_return_sum_mf_out((i-1)*n_x+j);
    end
end
% %/////////////   ���� �ֱ��� ////////////////////////
% hold on
% for i=1:1:n_y
%     
%     plot(abs(signal_return_sum_array(:,i)))
% end
% hold off

% figure(16)  %�������ϵ��Ӳ���
% hold on
% for i=1:1:n_y
%     
%     plot(abs(signal_return_1_array(:,i)+noise_return_1_array(:,i)))
% end
% hold off

% figure(7);mesh(1:n_y,1:n_x,real(signal_return_1_array));
% figure(6);mesh(1:n_y,1:n_x,abs(signal_return_1_array+  noise_return_1_array));
% figure(5);mesh(1:n_y,1:n_x,real(signal_return_sum_array));

%%%%%%%%%%%%%%%%%%%%%%%  ��FFT  �õ��ٶȷ���%%%%%%%%%%%%%%%%%%%%%
signal_return_1_array_FFT=zeros(n_x,n_y);
noise_return_1_array_FFT=zeros(n_x,n_y);
% signal_return_2_array_FFT=zeros(n_x,n_y);
signal_return_sum_array_FFT=zeros(n_x,n_y);
for i=1:n_x
    signal_return_1_array_FFT(i,:)=fft(signal_return_1_array(i,:),n_add);
    noise_return_1_array_FFT(i,:)=fft(noise_return_1_array(i,:),n_add);
    signal_return_sum_array_FFT(i,:)=fft(signal_return_sum_array(i,:),n_add);
end

%%%%%%%%%%%%%%%%%%%%%%%%% ����FFT������  %%%%%%%%%%%%%%%
[~, b_n]=max(signal_return_1_array(:,1));  %�õ���ѹ��������λ��
s_beforeFTT=signal_return_1_array(b_n,:);
n_beforeFTT=noise_return_1_array(b_n,:);
s_FFT_max=signal_return_1_array_FFT(b_n,:);  %�õ��ź�FFT����������
n_FFT_max=noise_return_1_array_FFT(b_n,:);   %�õ�����FFT����������

% [~, n_fft]=max(s_FFT_max);

%%%%fft��Ĵ���
signal_return_1_array_FFT_interp=interp(signal_return_1_array_FFT(b_n,:), 50);
f_fft=linspace(0,fs,length(signal_return_1_array_FFT_interp));
f_fft_2=linspace(0,PRF,length(signal_return_1_array_FFT_interp));
% plot(f_fft,10*log10(abs(signal_return_1_array_FFT_interp)),'linewidth',1);
% plot(f_fft_2,10*log10(abs(signal_return_1_array_FFT_interp)),'linewidth',1);

% �����ϵ�����Ϊ10*log10(N)  NΪFFT �ĵ���  ����Ϊ����ĸ��� n_add
[SNR_FFT_IN,~,~]=SNR_caj_time(s_beforeFTT,n_beforeFTT);
[SNR_FFT_out,~,~]=SNR_max_out(s_FFT_max,n_FFT_max);
Gain_fft=SNR_FFT_out-SNR_FFT_IN;

%%%%������ź��ٶ��ŵ�λ��
sort_A = sort(signal_return_sum_array(:,1), 'descend'); %��������
R_first_max = sort_A(1);
R_second_max = sort_A(2);
[~, R1_row] = find(signal_return_sum_array(:,1).' == R_first_max);    %������
[~, R2_row] = find(signal_return_sum_array(:,1).' == R_second_max);   

sort_B = sort(signal_return_sum_array_FFT(R1_row,:), 'descend'); %��������
v_first_max = sort_B(1);
v_second_max = sort_B(2);
[~, v1_row] = find(signal_return_sum_array_FFT(R1_row,:) == v_first_max);     %�ٶ���
[~, v2_row] = find(signal_return_sum_array_FFT(R1_row,:) == v_second_max );    

% %/////////////   ���� �ֱ��� ////////////////////////
% figure(17)
% hold on
% for i=1:1:n_y
%     
%     plot(abs(signal_return_sum_array_FFT(:,i))),xlim([2000,3000])
% end
% hold off

%/////////////  �ٶȷֱ��� ////////////////////////
figure(4)
hold on
for i=1:1:n_x
    
    plot(abs(signal_return_sum_array_FFT(i,:))),xlim([100,140])
end
hold off




% figure(15)  %�ٶ����ϵ��Ӳ���
% hold on
% for i=1:1:n_x
%     
%     plot(abs(signal_return_1_array_FFT(i,:)+noise_return_1_array_FFT(i,:)))
% end
% hold off



figure(1);mesh(1:n_y,1:n_x,abs(signal_return_1_array_FFT));
figure(2);mesh(1:n_y,1:n_x,real(signal_return_1_array_FFT+ noise_return_1_array_FFT));
figure(3);mesh(1:n_y,1:n_x,abs(signal_return_sum_array_FFT));



function  [Phase_code_signal_return]=Phase_code_signal_return_modulation(Phase_code_signal,fs,fc,Duty_Ratio,R_target,v_target,A_target,n_add)
    %�����״�ز��ź�
    %�Է����ź���ʱ���Ӷ�����Ƶ��
    c=3e8;    %����
    n_delay=round(2*R_target*fs/c);    %���ھ����������ʱ
    Phase_code_signal_delay=[zeros(1,n_delay),Phase_code_signal(1:length(Phase_code_signal)-n_delay)];
    fd=2*v_target*fc/c;           %���������Ƶ��
    t=1/fs:1/fs:length(Phase_code_signal_delay)/fs;  
    duppler=exp(1i*2*pi*fd*t);    %����������Ƶ������
    Phase_code_signal_return=A_target*Phase_code_signal_delay.*duppler;
end

 function [Phase_code_signal]=Phase_code_modulation(u,fd,fc,fs,Duty_Ratio,n_add)  
      %uΪm����   fdΪ��Ƶ    fcΪ��Ƶ   fsΪ����Ƶ��
      %����α�����λ�����ź�
   t = 0:1/fs:(1/fd*length(u)/Duty_Ratio*n_add-1/fs);      % һ����Ԫ��������ʱ���ڵĲ�����ʱ��
   code_signal= exp(1i*2*pi*fc*t);
   N=fs/fd;
   ck=[1, N*length(u)];
   for i=1:length(u)
       for j=1:N
           ck((i-1)*N+j)=u(i);
       end
   end
   ck_zeros=[ck,zeros(1,length(ck)*(1-Duty_Ratio)*10)];
    m_cpi= repmat(ck_zeros,1,n_add);
    Phase_code_signal=  m_cpi.*code_signal;
%    Phase_code_signal= ck_zeros.*code_signal;
%    Phase_code_signal=repmat(Phase_code_signal,1,n_add);

 end
 
 function [signal_demodulation]=demodulation(signal,fc,fs)
  %�Իز��ź� ���н��
   [~,b] = find(signal~=0);  
   n=length(signal);
   t=[zeros(1,b(1)-1),0:1/fs:(n-b)/fs];
   demodulate= exp(-1i*2*pi*fc*t);
   signal_demodulation=signal.*demodulate;
 end
 
 function [s,s_add_noise,SNR_in,xnoise]=add_SNR_xnoise(s_Preset,SNR_Preset,xnoise,s_length,N_filter)
    %  add_SNR_xnoise���ź��м���̶����������
    %  s_PresetΪ�������������ź�    SNR_PresetΪ�����  A_xnoiseΪ�����ı�׼��
    %  s_lengthΪ�������źŵĳ���
    %  sΪ�ź�     s_add_noise ������������ź�   SNR_in���������    xnoiseΪ���� 
    %  N_filterΪ�˲����Ľ���
    [~,b] = find(s_Preset~=0);    
%     xnoise = A_xnoise*randn(size(s_Preset))+1i*A_xnoise*randn(size(s_Preset));     %������ֵΪ0������ΪA_xnoiseƽ���ĸ�˹������
    Pn_add=sum(xnoise .*conj(xnoise ))/length(xnoise);                                % �������������ƽ������ 
    Ps_add=sum(s_Preset(b(1)+N_filter:b(1)+s_length-1).*conj(s_Preset(b(1)+N_filter:b(1)+s_length-1)))/(s_length-N_filter);            % ���㵥λ�źŵ�ƽ������ 
    A_1= sqrt((10^(SNR_Preset/10))*(Pn_add/Ps_add));     %�����Ӧ��������źŵķ���
    s=A_1.*s_Preset;                                    %�õ���Ӧ������µ��ź�
    s_add_noise=s+xnoise;                            % �����������ź�
    Ps_add_out=sum(s(b(1)+N_filter:b(1)+s_length-1).*conj(s(b(1)+N_filter:b(1)+s_length-1)))/(s_length-N_filter);
    SNR_in=10*log10(Ps_add_out/Pn_add);                   %���������
 end

  function [s,s_add_noise,SNR_in,xnoise]=add_SNR_xnoise_twotarget(s_Preset_1,s_Preset_2,SNR_Preset,A_xnoise,s_length)
    %  add_SNR_xnoise��˫�ز��ź��м���̶����������
    %  s_PresetΪ�������������ź�    SNR_PresetΪ�����  A_xnoiseΪ�����ı�׼��
    %  s_lengthΪ�������źŵĳ���
    %  sΪ�ź�     s_add_noise ������������ź�   SNR_in���������    xnoiseΪ���� 
    [~,b] = find(s_Preset_1~=0);    
    xnoise_1 = A_xnoise*randn(size(s_Preset_1))+ 1i*A_xnoise*randn(size(s_Preset_1));     %������ֵΪ0������ΪA_xnoiseƽ���ĸ�˹������
    xnoise_2 = A_xnoise*randn(size(s_Preset_2))+ 1i*A_xnoise*randn(size(s_Preset_1));     %������ֵΪ0������ΪA_xnoiseƽ���ĸ�˹������
    if  length(xnoise_1 )>=length(xnoise_2 )
        xnoise=xnoise_1;
    else
        xnoise=xnoise_2;
    end
    Pn_add=sum(xnoise .*conj(xnoise ))/length(xnoise);                                % �������������ƽ������ 
    Ps_add_1=sum(s_Preset_1(b(1):b(1)+s_length-1).*conj(s_Preset_1(b(1):b(1)+s_length-1)))/(s_length);            % ���㵥λ�źŵ�ƽ������ 
    [~,a] = find(s_Preset_2~=0); 
    Ps_add_2=sum(s_Preset_2(a(1):a(1)+s_length-1).*conj(s_Preset_2(a(1):a(1)+s_length-1)))/(s_length);            % ���㵥λ�źŵ�ƽ������ 
    Ps_add=(Ps_add_1+Ps_add_2)/2;
    A_1= sqrt((10^(SNR_Preset/10))*(Pn_add/Ps_add));     %�����Ӧ��������źŵķ���
     %ʵ��������ͬ���ȵ��ź����
      r1=length(s_Preset_1);
      r2=length(s_Preset_2);
      s_sum=zeros(1,max(r1,r2));
      s_sum(1,1:r1)=s_Preset_1;
     
      %�ϳɻز��ź�
      s_sum(1,1:r2)=s_sum(1,1:r2)+s_Preset_2;
      s=A_1*s_sum;
      s_add_noise=s+xnoise;                            % �����������ź�
      Ps_add_out=A_1^2*(Ps_add_1+Ps_add_2)/2;
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
   SK_out=SK.*HK;                %LFM�ź���MF��FFTȻ�����
   s_mf_out=ifft(SK_out,N_s);    %�õ���ѹ����ź�
  
   N_n=2^ceil(log2(N3+N2-1));
   NK=fft(xnoise_in,N_n);
   HK_n=fft(h_in,N_n);
   NK_out=NK.*HK_n; 
   xnoise_mf_out_temp=ifft(NK_out,N_n);    %�õ���ѹ�������
   xnoise_mf_out=xnoise_mf_out_temp(N_n/2-N_n/2*(N_s/N_n)+1:N_n/2+N_n/2*(N_s/N_n));  %�������нض�
  end

function  [SNR_out,p_s_max,p_n]=SNR_max_out(s,xnoise)
  %���������ȼ��źŷ�ֵ���ʡ�����ƽ������
  p_s_max=10*log10(max(s)*conj(max(s)));
  p_n=10*log10(sum(xnoise.*conj(xnoise))/length(xnoise));
  SNR_out =p_s_max- p_n;
end



function  [SNR_out,p_s,p_n]=SNR_caj_time(s,xnoise)
  %������ȼ�����
  p_s=10*log10(sum(s.*conj(s))/length(s));
  p_n=10*log10(sum(xnoise.*conj(xnoise))/length(xnoise));
  SNR_out =p_s- p_n;
end

function  [signal_out]=filter_compose(signal,A_filter,b_filter,s_length,Duty_Ratio,N_add)
    [~,b] = find(signal~=0); 
      signal_temp=zeros(N_add,s_length);
     for i=1:N_add
      signal_temp(i,:)=[filter(b_filter,A_filter,signal(b(1)+(i-1)*s_length:(i-1)*s_length+b(1)+s_length*Duty_Ratio-1)),zeros(1,s_length*(1-Duty_Ratio))];
     end
     
      signal_out=zeros(1,b(1)-1);
     for i=1:N_add
         signal_out=[signal_out,signal_temp(i,:)];
     end
 
%      signal_out=[zeros(1,b-1), signal_out];
     signal_out=signal_out(1:length(signal));
end

function [s,s_add_noise,SNR_in,xnoise]=add_SNR(s_Preset,SNR_Preset,xnoise)
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


