clear
% ��ȡ���м����źţ�������ά����ÿ��Ϊһ���ļ��е���������*����ͨ��
% Ԥ�������
data=zeros(15,27000,3);
% ��ȡ
for i=1:15
    data((i),:,:)=load_emg_from_txt_and_remove_relax_1_5(['D:\nutstore\TCDS_Special_Issue\Code_and_Data\emg_txt\kf\' num2str(i) '.txt'],1);
end
% ���䱣��������Ŀ¼
% save('emg_data_kf.mat','data')