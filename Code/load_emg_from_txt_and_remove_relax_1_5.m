%������emg�źŶ�ȡ���������ÿ��Ի��ʱ�䡢ͨ��1��ͨ��2
%ʹ�÷��� data=load_emg('emg/1-1.txt',1);
%�����Ҫ��ͼ����Ҫ����1���������Ҫ��Ϊ0
function data=load_emg_from_txt_and_remove_relax_1_5(filename,plot_flag)
%% ��ʼ��
%countΪ������������ѭ���и�������
count=1;
%�����վ���
data=zeros(27000,3);
raw_data=load (filename);
%nΪ���ݵĳ���
n=max(size(raw_data));
t=raw_data(:,1);
channel1=raw_data(:,2);
channel2=raw_data(:,3);
% channel3=raw_data(:,4);
%% ��ȡ
%ѭ�������������źŵ㣬�����е������źŴ���data����
%data����Ϊ3�еı���������ÿ��������3s*1500Hz���źŵ�
%����3s*1500Hz*6����=27000�У�3��
%��ʱ����ȡ��ȥ��ÿ���źŵ�ͷβ1s
for i=1:n
    if(t(i)>6&&t(i)<=9)
       data(count,1)=t(i);
       data(count,2)=channel1(i);
       data(count,3)=channel2(i);
       count=count+1;
    end
    if(t(i)>16&&t(i)<=19)
       data(count,1)=t(i);
       data(count,2)=channel1(i);
       data(count,3)=channel2(i);
       count=count+1;
    end
    if(t(i)>26&&t(i)<=29)
       data(count,1)=t(i);
       data(count,2)=channel1(i);
       data(count,3)=channel2(i);
       count=count+1; 
    end            
    if(t(i)>36&&t(i)<=39)
       data(count,1)=t(i);
       data(count,2)=channel1(i);
       data(count,3)=channel2(i);
       count=count+1;
    end
    if(t(i)>46&&t(i)<=49)
       data(count,1)=t(i);
       data(count,2)=channel1(i);
       data(count,3)=channel2(i);
       count=count+1;
    end
    if(t(i)>56&&t(i)<=59)
       data(count,1)=t(i);
       data(count,2)=channel1(i);
       data(count,3)=channel2(i);
       count=count+1;
    end
end
%% ��ͼ
if(plot_flag)
    figure('position',[300 200 800 500])
    for i=1:6
        for j=1:6
            if(rem(i,2))
                if(rem(j,2))
                    plot_color=[0 0 1];
                else
                    plot_color=[0 0 1];
                end
            else
                if(rem(j,2))
                    plot_color=[0 0 1];
                else
                    plot_color=[0 0 1];
                end
            end
            subplot(2,1,1)
            hold on
            plot(data(((i-1)*4500+(j-1)*750+1):((i-1)*4500+(j)*750),1),data(((i-1)*4500+(j-1)*750+1):((i-1)*4500+(j)*750),2),'Color',plot_color)
            axis([0;60;-1500;1500])
            subplot(2,1,2)
            hold on
            plot(data(((i-1)*4500+(j-1)*750+1):((i-1)*4500+(j)*750),1),data(((i-1)*4500+(j-1)*750+1):((i-1)*4500+(j)*750),3),'Color',plot_color)
            axis([0;60;-1500;1500])
        end
    end
end
% % ��ͼ
% if(plot_flag)
%     figure
%     for i=1:6
%         for j=1:6
%             if(rem(i,2))
%                 if(rem(j,2))
%                     plot_color=[1 0 0];
%                 else
%                     plot_color=[1 0 1];
%                 end
%             else
%                 if(rem(j,2))
%                     plot_color=[0 0 1];
%                 else
%                     plot_color=[0 0.5 1];
%                 end
%             end
%             subplot(2,1,1)
%             hold on
%             plot(((i-1)*4500+(j-1)*750+1):((i-1)*4500+(j)*750),...
%                 data(((i-1)*4500+(j-1)*750+1):((i-1)*4500+(j)*750),2),'Color',plot_color)
%             subplot(2,1,2)
%             hold on
%             plot(((i-1)*4500+(j-1)*750+1):((i-1)*4500+(j)*750),...
%                 data(((i-1)*4500+(j-1)*750+1):((i-1)*4500+(j)*750),3),'Color',plot_color)
%         end
%     end
% end
end