clear
load('kf_nosliding.mat')
%load('D:\nutstore\TCDS_Special_Issue\Code_and_Data\emg_nosliding\kf_nosliding.mat')
for experiment_=1:15%第几次实验
    for action_=1:6%第几个动作
        for window_=1:141%第几个窗口
            data(experiment_,action_,window_,:,1)=emg(experiment_,(action_-1)*4500+(window_-1)*30+1:(action_-1)*4500+(window_-1)*30+300,1);
            data(experiment_,action_,window_,:,2)=emg(experiment_,(action_-1)*4500+(window_-1)*30+1:(action_-1)*4500+(window_-1)*30+300,2);
        end
    end
end
save('kf_sliding.mat','data')
% figure
% plot(emg(5,4501:5000,1))
% axis([0;500;-100;100])
% figure
% for i=1:5
%     subplot(5,1,i)
%     d=reshape(data(5,2,i,:,1),[300,1]);
%     plot(1+(i-1)*30:(i-1)*30+300,d)
%     axis([0;500;-100;100])
% end