clear,clc
%accuracy_ada,accuracy_bp:����*����+ƽ��*1ʱ��2Ƶ��3ʱ��+Ƶ��/��ȷ��
%fd_ada,fd_bp,td_ada,td_bp,tfd_ada,tfd_bp:����*635*����/ʶ����
%train_time_ada,train_time_bp:����*1ʱ��2Ƶ��3ʱ��+Ƶ��/ѵ��ʱ��
%test_time_ada,test_time_bp:����*1ʱ��2Ƶ��3ʱ��+Ƶ��/����ʱ��
testname='yyk'
load(['D:\nutstore\TCDS_Special_Issue\Code_and_Data\half-results\' testname '_accuracy_ada.mat'])
load(['D:\nutstore\TCDS_Special_Issue\Code_and_Data\half-results\' testname '_accuracy_bp.mat'])
load(['D:\nutstore\TCDS_Special_Issue\Code_and_Data\half-results\' testname '_fd_ada.mat'])
load(['D:\nutstore\TCDS_Special_Issue\Code_and_Data\half-results\' testname '_fd_bp.mat'])
load(['D:\nutstore\TCDS_Special_Issue\Code_and_Data\half-results\' testname '_td_ada.mat'])
load(['D:\nutstore\TCDS_Special_Issue\Code_and_Data\half-results\' testname '_td_bp.mat'])
load(['D:\nutstore\TCDS_Special_Issue\Code_and_Data\half-results\' testname '_tfd_ada.mat'])
load(['D:\nutstore\TCDS_Special_Issue\Code_and_Data\half-results\' testname '_tfd_bp.mat'])
load(['D:\nutstore\TCDS_Special_Issue\Code_and_Data\half-results\' testname '_train_time_ada.mat'])
load(['D:\nutstore\TCDS_Special_Issue\Code_and_Data\half-results\' testname '_train_time_bp.mat'])
load(['D:\nutstore\TCDS_Special_Issue\Code_and_Data\half-results\' testname '_test_time_ada.mat'])
load(['D:\nutstore\TCDS_Special_Issue\Code_and_Data\half-results\' testname '_test_time_bp.mat'])
matrix_fd_ada=zeros(6,6);
matrix_fd_bp=zeros(6,6);
matrix_td_ada=zeros(6,6);
matrix_td_bp=zeros(6,6);
matrix_tfd_ada=zeros(6,6);
matrix_tfd_bp=zeros(6,6);
%ƽ����ȷ��
accu_ada=mean(accuracy_ada);
accu_bp=mean(accuracy_bp);
%ƽ��ʱ��
tr_time_ada=mean(train_time_ada)
tr_time_bp=mean(train_time_bp)
te_time_ada=mean(test_time_ada)
te_time_bp=mean(test_time_bp)
%��������
for i=1:100%����
     for j=1:6%����
          for k=1:635%ʶ�����ĸ���
               if(td_bp(i,k,j)>0.5&&td_bp(i,k,j)<1.5)
                    matrix_td_bp(j,1)=matrix_td_bp(j,1)+1;
               end
               if(td_bp(i,k,j)>1.5&&td_bp(i,k,j)<2.5)
                    matrix_td_bp(j,2)=matrix_td_bp(j,2)+1;
               end
               if(td_bp(i,k,j)>2.5&&td_bp(i,k,j)<3.5)
                    matrix_td_bp(j,3)=matrix_td_bp(j,3)+1;
               end
               if(td_bp(i,k,j)>3.5&&td_bp(i,k,j)<4.5)
                    matrix_td_bp(j,4)=matrix_td_bp(j,4)+1;
               end
               if(td_bp(i,k,j)>4.5&&td_bp(i,k,j)<5.5)
                    matrix_td_bp(j,5)=matrix_td_bp(j,5)+1;
               end
               if(td_bp(i,k,j)>5.5&&td_bp(i,k,j)<6.5)
                    matrix_td_bp(j,6)=matrix_td_bp(j,6)+1;
               end
          end
     end
end
