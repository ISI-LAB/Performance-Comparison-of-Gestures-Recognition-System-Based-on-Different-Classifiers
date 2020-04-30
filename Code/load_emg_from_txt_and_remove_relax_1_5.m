%基本的emg信号读取函数，调用可以获得时间、通道1和通道2
%使用方法 data=load_emg('emg/1-1.txt',1);
%如果需要画图，需要参数1，如果不需要则为0
function data=load_emg_from_txt_and_remove_relax_1_5(filename,plot_flag)
%% 初始化
%count为辅助变量，在循环中辅助计数
count=1;
%创建空矩阵
data=zeros(27000,3);
raw_data=load (filename);
%n为数据的长度
n=max(size(raw_data));
t=raw_data(:,1);
channel1=raw_data(:,2);
channel2=raw_data(:,3);
% channel3=raw_data(:,4);
%% 提取
%循环，遍历所有信号点，将其中的有用信号存入data变量
%data变量为3列的变量，其中每个动作有3s*1500Hz个信号点
%共计3s*1500Hz*6动作=27000行，3列
%按时间提取，去掉每个信号的头尾1s
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
%% 画图
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
% % 画图
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