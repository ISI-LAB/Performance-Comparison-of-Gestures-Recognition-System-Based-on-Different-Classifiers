clear
% 读取所有肌电信号，生成三维矩阵。每层为一个文件中的六个动作*两个通道
% 预分配变量
data=zeros(15,27000,3);
% 读取
for i=1:15
    data((i),:,:)=load_emg_from_txt_and_remove_relax_1_5(['D:\nutstore\TCDS_Special_Issue\Code_and_Data\emg_txt\kf\' num2str(i) '.txt'],1);
end
% 将其保存至工作目录
% save('emg_data_kf.mat','data')