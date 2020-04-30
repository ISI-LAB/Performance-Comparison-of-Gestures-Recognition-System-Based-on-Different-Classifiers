%本程序可以提取特征，并且计算每种特征所用的提取时间
clear
clc
testname='lbh'
load(['D:\nutstore\TCDS_Special_Issue\Code_and_Data\emg_sliding\' testname '_sliding.mat'])
%%时域特征
rms_feature=zeros(15,6,141,2);
rms_time=zeros(100,1);
for count_=1:100%计算次数
    tic%计算rms时间
    for experiment_=1:15%对于每次实验
        for action_=1:6%对于每个动作
            for window_=1:141%对于每个窗口
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
for count_=1:100%计算次数
    tic%计算zc时间
    for experiment_=1:15%对于每次实验
        for action_=1:6%对于每个动作
            for window_=1:141%对于每个窗口
                for i=1:300-1%对于每个数据点
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
for count_=1:100%计算次数
    tic%计算wamp时间
    for experiment_=1:15%对于每次实验
        for action_=1:6%对于每个动作
            for window_=1:141%对于每个窗口
                for i=1:300-1%对于每个数据点
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

%%频域特征MPS
mps_feature=zeros(15,6,141,2);
mps_time=zeros(100,1);
Fs = 1500;			% 信号采样频率
T = 1/Fs;			% 信号采样周期
L = 300;			% 信号长度
t = (0:L-1)*T;		% 时间向量
f = Fs*(0:(L/2))/L;
Iv = 1:length(f);	% Index Vector
for count_=1:100%计算次数
	tic
    for experiment_=1:15%对于每次实验
        for action_=1:6%对于每个动作
            for window_=1:141%对于每个窗口
            	%傅里叶变换开始
            	%波形 S1=通道1（尺侧腕伸肌） S2=通道2（桡侧腕屈肌）
		        S1 = reshape(data(experiment_,action_,window_,:,1),[1 300]);
		        S2 = reshape(data(experiment_,action_,window_,:,2),[1 300]);
		        %变换
		        Y1 = fft(S1);
		        Y2 = fft(S2);
		        %对称变单侧，乘2
		        P1_2 = abs(Y1/L);
		        P1_1 = P1_2(1:fix(L/2)+1);
		        P2_2 = abs(Y2/L);
		        P2_1 = P2_2(1:fix(L/2)+1);
		        %P1表示频谱值
		        P1_1(2:end-1) = 2*P1_1(2:end-1);
		        P2_1(2:end-1) = 2*P2_1(2:end-1);
		        %傅里叶变换完毕
		        %提取MPS     mi=maxindex
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
%%频域特征MF
mf_feature=zeros(15,6,141,2);
mf_time=zeros(100,1);
Fs = 1500;			% 信号采样频率
T = 1/Fs;			% 信号采样周期
L = 300;			% 信号长度
t = (0:L-1)*T;		% 时间向量
f = Fs*(0:(L/2))/L;
Iv = 1:length(f);	% Index Vector
for count_=1:100%计算次数
	tic
    for experiment_=1:15%对于每次实验
        for action_=1:6%对于每个动作
            for window_=1:141%对于每个窗口
            	%傅里叶变换开始
            	%波形 S1=通道1（尺侧腕伸肌） S2=通道2（桡侧腕屈肌）
		        S1 = reshape(data(experiment_,action_,window_,:,1),[1 300]);
		        S2 = reshape(data(experiment_,action_,window_,:,2),[1 300]);
		        %变换
		        Y1 = fft(S1);
		        Y2 = fft(S2);
		        %对称变单侧，乘2
		        P1_2 = abs(Y1/L);
		        P1_1 = P1_2(1:fix(L/2)+1);
		        P2_2 = abs(Y2/L);
		        P2_1 = P2_2(1:fix(L/2)+1);
		        %P1表示频谱值
		        P1_1(2:end-1) = 2*P1_1(2:end-1);
		        P2_1(2:end-1) = 2*P2_1(2:end-1);
		        %傅里叶变换完毕
		        %提取MF
		        CumAmp1 = cumtrapz(f, abs(P1_1(Iv))); % Integrate FFT Amplitude
		        mf_feature(experiment_,action_,window_,1)=interp1(CumAmp1, f, CumAmp1(end)/2);% Use‘interp1’To Find‘MF’
		        CumAmp2 = cumtrapz(f, abs(P2_1(Iv))); % Integrate FFT Amplitude
		        mf_feature(experiment_,action_,window_,2)=interp1(CumAmp2, f, CumAmp2(end)/2);% Use‘interp1’To Find‘MF’
            end
        end
    end
	mf_time(count_)=toc;
    disp(['mf' num2str(count_)])
end
disp('mf_done!')
mpf_feature=zeros(15,6,141,2);
mpf_time=zeros(100,1);
Fs = 1500;			% 信号采样频率
T = 1/Fs;			% 信号采样周期
L = 300;			% 信号长度
t = (0:L-1)*T;		% 时间向量
f = Fs*(0:(L/2))/L;
Iv = 1:length(f);	% Index Vector
for count_=1:100%计算次数
	tic
    for experiment_=1:15%对于每次实验
        for action_=1:6%对于每个动作
            for window_=1:141%对于每个窗口
            	%傅里叶变换开始
            	%波形 S1=通道1（尺侧腕伸肌） S2=通道2（桡侧腕屈肌）
		        S1 = reshape(data(experiment_,action_,window_,:,1),[1 300]);
		        S2 = reshape(data(experiment_,action_,window_,:,2),[1 300]);
		        %提取MPF
		        %[pxx,f] = pwelch(x,window,noverlap,f,fs)
		        %x=输入信号 window=窗口大小 可以为空 noverlap=重叠数 可以为空 f=周期频率 fs=采样率
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