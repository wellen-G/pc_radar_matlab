%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
% ������������������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc 
% clear 
close all 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  m_length=127;                     %�볤Ϊ127
fm=72e6;                          %��ƵΪ72MHz
fc=10e9;                          %�ź���ƵΪ 10GHz
fs=72e6*4;                         %������Ϊ200MHz
Ts=1/fs;
Duty_Ratio= 0.10;                  %ռ�ձ�Ϊ10%
PRT=m_length*(1/fm)/Duty_Ratio;    %�����ظ�����
PRF=1/PRT;                        %�����ظ�Ƶ��
n_add=round(100e-3/(m_length/fm)*Duty_Ratio);        %����ۼƸ���
n_add=128;
CPI=n_add*PRT;                       %��δ���ʱ��
SNR_Preset=-10;                       %���������
A_xnoise =1;                            %������׼��
N_PRT=m_length*(fs/fm)/Duty_Ratio;    %һ�������ڵĲ�������

T_width=m_length/fm;              %�źŵ�ʱ��
B=fm;
%%%%%%%%%%%%%%%%%%%%%%  ����m����   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ��˹����ź� 
u = idinput(m_length, 'prbs',[0,1],[-1,1])';   %����m����
%%%%%%%%%%%%%%%%%%%%%%%% ����α�����λ��������ź� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
Phase_code_signal=Phase_code_modulation( u,fm,fc,fs,Duty_Ratio,n_add) ;      %���������ź�
%    plot( real(Phase_code_signal)),   
s_length=length(Phase_code_signal);                        %�����źŵĳ���

%%%%%%%%%%%%%%%%%%%%%%%%%%%  �����ز��ź� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R_target=[0,14200];               %Ŀ��ľ���    R_max=(3e8/2)*PRT*(1-Duty_Ratio)=14288  R_min=(3e8/2)*PRT*Duty_Ratio=1587
v_target=[30,100];                    %Ŀ����ٶ�    c_min=3e8*PRF/2/fc=141.7
A_target=[1,0.3];                  %Ŀ��ز��źŵķ�ֵ˥��
fd=linspace(0,1/(127/fm)*5,300);
mf_max=zeros(1,length(fd));          %��ѹ��������ֵ
mf_max_2=zeros(1,length(fd));        %��ѹ��������ֵ
mf_sec=zeros(1,length(fd));          %��ѹ����԰��ֵ
max_sec=zeros(1,length(fd));         %���԰��ֵ
b_max_zhuban=zeros(1,length(fd));  %������ѹ�����ֵ��λ�á�
for i=1:length(fd)

   Phase_code_signal_return_1=Phase_code_signal_return_modulation(Phase_code_signal,fs,fc,Duty_Ratio,R_target(1),v_target(1),A_target(1),fd(i));
   h=conj(fliplr(Phase_code_signal(1:m_length*round(fs/fm))));   %����ƥ���˲������Գ弤��Ӧ 
   mf_out=conv(Phase_code_signal_return_1(1:m_length*round(fs/fm)),h);
%    [SLL,S_sec]=SLL_CAL(mf_out);
%   mf_sec(i)=S_sec;
   [mf_max(i),b_max_zhuban(i)]=max(abs(mf_out));
   mf_max_2(i)=abs(mf_out(b_max_zhuban(1)));
%    max_sec(i)=mf_max(i)/mf_sec(i);
end

%  plot(fd,20*log10(abs(mf_max)),'linewidth',2),title('��������������')
%  plot(fd,20*log10(abs(mf_max_2)),'linewidth',2),title('��������������')
%  %��������λ����ͬ�µĶ��������ޡ�



function  [Phase_code_signal_return]=Phase_code_signal_return_modulation(Phase_code_signal,fs,fc,Duty_Ratio,R_target,v_target,A_target,fd)
    %�����״�ز��ź�
    %�Է����ź���ʱ���Ӷ�����Ƶ��
    c=3e8;    %����
    n_delay=round(2*R_target*fs/c);    %���ھ����������ʱ
    Phase_code_signal_delay=[zeros(1,n_delay),Phase_code_signal(1:length(Phase_code_signal)-n_delay)];
%     fd=2*v_target*fc/c;           %���������Ƶ��
    t=1/fs:1/fs:length(Phase_code_signal_delay)/fs;  
    duppler=exp(1i*2*pi*fd*t);    %����������Ƶ������
    Phase_code_signal_return=A_target*Phase_code_signal_delay.*duppler;
end

 function [Phase_code_signal]=Phase_code_modulation(u,fd,fc,fs,Duty_Ratio,n_add)     %����һ�������ڵ�α�����λ�����ź�
   %uΪm����   fdΪ��Ƶ    fcΪ��Ƶ   fsΪ����Ƶ��
   
   t = 0:1/fs:(1/fd*length(u)/Duty_Ratio-1/fs);      % һ����Ԫ��������ʱ���ڵĲ�����ʱ��
   code_signal= exp(1i*2*pi*fc*t);
   N=fs/fd;
   ck=[1, N*length(u)];
   for i=1:length(u)
       for j=1:N
           ck((i-1)*N+j)=u(i);
       end
   end
   ck_zeros=[ck,zeros(1,length(ck)*(1-Duty_Ratio)*10)];
%    m_cpi= repmat(ck_zeos,1,n_add);
   Phase_code_signal= ck_zeros;%.*code_signal;  %������Ƶ
   Phase_code_signal=repmat(Phase_code_signal,1,n_add);

 end
 
 function [SLL,S_sec]=SLL_CAL(signal)
    %���źŵ��԰��ƽ
   
%     signal=signal(round(11/30*length(signal)):round(19/30*length(signal)));
    Lmax=diff(sign(diff(abs(signal))))==-2;    %�õ���ֵ��
    Lmax=[false,Lmax,false];
    vmax=signal(Lmax);
    [S_max,S_max_index]=max(vmax);             %��������ƽ��λ��
    S_sec_sum=[vmax(S_max_index-3:S_max_index-1),vmax(S_max_index+1:S_max_index+3)];  %�õ�������Χ���԰�
    S_sec=max(S_sec_sum);               %�õ�����԰�
%     Signal_sort=sort(abs(vmax));    %����������
%     S_max= Signal_sort(length(Signal_sort));
%     S_sec=Signal_sort(length(Signal_sort)-1);
    SLL=20*log10(abs(S_sec)/abs(S_max));
 end