%�����������ȡ���������Ҽ���ÿ���������õ���ȡʱ��
clear
clc
testname='lbh'
load(['D:\nutstore\TCDS_Special_Issue\Code_and_Data\emg_sliding\' testname '_sliding.mat'])
%%ʱ������
rms_feature=zeros(15,6,141,2);
rms_time=zeros(100,1);
for count_=1:100%�������
    tic%����rmsʱ��
    for experiment_=1:15%����ÿ��ʵ��
        for action_=1:6%����ÿ������
            for window_=1:141%����ÿ������
                rms_feature(experiment_,action_,window_,1)=log10(rms(data(experiment_,action_,window_,:,1)));
                rms_feature(experiment_,action_,window_,2)=log10(rms(data(experiment_,action_,window_,:,2)));
            end
        end
    end
    rms_time(count_)=toc;
    if rem(count_,10)==0
        disp(['rms' num2str(count_)])
    end
end
disp('rms_done!')
zc_feature=zeros(15,6,141,2);
zc_time=zeros(100,1);
for count_=1:100%�������
    tic%����zcʱ��
    for experiment_=1:15%����ÿ��ʵ��
        for action_=1:6%����ÿ������
            for window_=1:141%����ÿ������
                for i=1:300-1%����ÿ�����ݵ�
                    if ((data(experiment_,action_,window_,i,1)*data(experiment_,action_,window_,i+1,1))<0)
                        zc_feature(experiment_,action_,window_,1)=zc_feature(experiment_,action_,window_,1)+1;
                    end
                    if ((data(experiment_,action_,window_,i,2)*data(experiment_,action_,window_,i+1,2))<0)
                        zc_feature(experiment_,action_,window_,2)=zc_feature(experiment_,action_,window_,2)+1;
                    end
                end
            end
        end
    end
    zc_time(count_)=toc;
    if rem(count_,10)==0
        disp(['zc' num2str(count_)])
    end
end
disp('zc_done!')
wamp_feature=zeros(15,6,141,2);
wamp_time=zeros(100,1);
wamp_threshold=100;
for count_=1:100%�������
    tic%����wampʱ��
    for experiment_=1:15%����ÿ��ʵ��
        for action_=1:6%����ÿ������
            for window_=1:141%����ÿ������
                for i=1:300-1%����ÿ�����ݵ�
                    if (abs((data(experiment_,action_,window_,i,1)-data(experiment_,action_,window_,i+1,1)))>wamp_threshold)
                        wamp_feature(experiment_,action_,window_,1)=wamp_feature(experiment_,action_,window_,1)+1;
                    end
                    if (abs((data(experiment_,action_,window_,i,2)-data(experiment_,action_,window_,i+1,2)))>wamp_threshold)
                        wamp_feature(experiment_,action_,window_,2)=wamp_feature(experiment_,action_,window_,2)+1;
                    end
                end
            end
        end
    end
    wamp_time(count_)=toc;
    if rem(count_,10)==0
        disp(['wamp' num2str(count_)])
    end
end
disp('wamp_done!')

%%Ƶ������MPS
mps_feature=zeros(15,6,141,2);
mps_time=zeros(100,1);
Fs = 1500;			% �źŲ���Ƶ��
T = 1/Fs;			% �źŲ�������
L = 300;			% �źų���
t = (0:L-1)*T;		% ʱ������
f = Fs*(0:(L/2))/L;
Iv = 1:length(f);	% Index Vector
for count_=1:100%�������
	tic
    for experiment_=1:15%����ÿ��ʵ��
        for action_=1:6%����ÿ������
            for window_=1:141%����ÿ������
            	%����Ҷ�任��ʼ
            	%���� S1=ͨ��1���߲����켡�� S2=ͨ��2�������������
		        S1 = reshape(data(experiment_,action_,window_,:,1),[1 300]);
		        S2 = reshape(data(experiment_,action_,window_,:,2),[1 300]);
		        %�任
		        Y1 = fft(S1);
		        Y2 = fft(S2);
		        %�ԳƱ䵥�࣬��2
		        P1_2 = abs(Y1/L);
		        P1_1 = P1_2(1:fix(L/2)+1);
		        P2_2 = abs(Y2/L);
		        P2_1 = P2_2(1:fix(L/2)+1);
		        %P1��ʾƵ��ֵ
		        P1_1(2:end-1) = 2*P1_1(2:end-1);
		        P2_1(2:end-1) = 2*P2_1(2:end-1);
		        %����Ҷ�任���
		        %��ȡMPS     mi=maxindex
		        mps_feature(experiment_,action_,window_,1)=max(P1_1);
		        mps_feature(experiment_,action_,window_,2)=max(P2_1);
            end
        end
    end
	mps_time(count_)=toc;
    if rem(count_,10)==0
        disp(['mps' num2str(count_)])
    end
end
disp('mps_done!')
%%Ƶ������MF
mf_feature=zeros(15,6,141,2);
mf_time=zeros(100,1);
Fs = 1500;			% �źŲ���Ƶ��
T = 1/Fs;			% �źŲ�������
L = 300;			% �źų���
t = (0:L-1)*T;		% ʱ������
f = Fs*(0:(L/2))/L;
Iv = 1:length(f);	% Index Vector
for count_=1:100%�������
	tic
    for experiment_=1:15%����ÿ��ʵ��
        for action_=1:6%����ÿ������
            for window_=1:141%����ÿ������
            	%����Ҷ�任��ʼ
            	%���� S1=ͨ��1���߲����켡�� S2=ͨ��2�������������
		        S1 = reshape(data(experiment_,action_,window_,:,1),[1 300]);
		        S2 = reshape(data(experiment_,action_,window_,:,2),[1 300]);
		        %�任
		        Y1 = fft(S1);
		        Y2 = fft(S2);
		        %�ԳƱ䵥�࣬��2
		        P1_2 = abs(Y1/L);
		        P1_1 = P1_2(1:fix(L/2)+1);
		        P2_2 = abs(Y2/L);
		        P2_1 = P2_2(1:fix(L/2)+1);
		        %P1��ʾƵ��ֵ
		        P1_1(2:end-1) = 2*P1_1(2:end-1);
		        P2_1(2:end-1) = 2*P2_1(2:end-1);
		        %����Ҷ�任���
		        %��ȡMF
		        CumAmp1 = cumtrapz(f, abs(P1_1(Iv))); % Integrate FFT Amplitude
		        mf_feature(experiment_,action_,window_,1)=interp1(CumAmp1, f, CumAmp1(end)/2);% Use��interp1��To Find��MF��
		        CumAmp2 = cumtrapz(f, abs(P2_1(Iv))); % Integrate FFT Amplitude
		        mf_feature(experiment_,action_,window_,2)=interp1(CumAmp2, f, CumAmp2(end)/2);% Use��interp1��To Find��MF��
            end
        end
    end
	mf_time(count_)=toc;
    disp(['mf' num2str(count_)])
end
disp('mf_done!')
mpf_feature=zeros(15,6,141,2);
mpf_time=zeros(100,1);
Fs = 1500;			% �źŲ���Ƶ��
T = 1/Fs;			% �źŲ�������
L = 300;			% �źų���
t = (0:L-1)*T;		% ʱ������
f = Fs*(0:(L/2))/L;
Iv = 1:length(f);	% Index Vector
for count_=1:100%�������
	tic
    for experiment_=1:15%����ÿ��ʵ��
        for action_=1:6%����ÿ������
            for window_=1:141%����ÿ������
            	%����Ҷ�任��ʼ
            	%���� S1=ͨ��1���߲����켡�� S2=ͨ��2�������������
		        S1 = reshape(data(experiment_,action_,window_,:,1),[1 300]);
		        S2 = reshape(data(experiment_,action_,window_,:,2),[1 300]);
		        %��ȡMPF
		        %[pxx,f] = pwelch(x,window,noverlap,f,fs)
		        %x=�����ź� window=���ڴ�С ����Ϊ�� noverlap=�ص��� ����Ϊ�� f=����Ƶ�� fs=������
		        [pxx1,f_1] = pwelch(S1,[],[],500,Fs);
		        mpf_feature(experiment_,action_,window_,1)=sum(f_1.*pxx1)/sum(pxx1);
		        [pxx2,f_2] = pwelch(S2,[],[],500,Fs);
		        mpf_feature(experiment_,action_,window_,2)=sum(f_2.*pxx2)/sum(pxx2);
            end
        end
    end
	mpf_time(count_)=toc;
    disp(['mpf' num2str(count_)])
end
disp('mpf_done!')

save([testname '_zc.mat'],'zc_feature')
save([testname '_rms.mat'],'rms_feature')
save([testname '_wamp.mat'],'wamp_feature')
save([testname '_mpf.mat'],'mpf_feature')
save([testname '_mf.mat'],'mf_feature')
save([testname '_mps.mat'],'mps_feature')
save([testname '_zc_time.mat'],'zc_time')
save([testname '_rms_time.mat'],'rms_time')
save([testname '_wamp_time.mat'],'wamp_time')
save([testname '_mpf_time.mat'],'mpf_time')
save([testname '_mf_time.mat'],'mf_time')
save([testname '_mps_time.mat'],'mps_time')

disp('all_done!')