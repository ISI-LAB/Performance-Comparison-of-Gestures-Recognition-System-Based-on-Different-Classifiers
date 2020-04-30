clear
load('emg_data_kf')
for experiment_=1:15%第几次实验
    emg(experiment_,:,1)=data(experiment_,:,2);
    emg(experiment_,:,2)=data(experiment_,:,3);
end
save('kf_nosliding.mat','emg')