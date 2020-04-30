%这个程序可以计算出BP神经网络和Adaboost的分类结果、分类所用时间
clear;clc;close all;warning off;
%样本姓名
testname='yyk';
%calculition_time循环次数
c_t=1;
%读取数据是按照分离度的从小到大排序的
load (['D:\OneDrive - mail.nankai.edu.cn\论文\201907TCDS\Code_and_Data\code\' testname '_rms.mat'])
load (['D:\OneDrive - mail.nankai.edu.cn\论文\201907TCDS\Code_and_Data\code\' testname '_zc.mat'])
load (['D:\OneDrive - mail.nankai.edu.cn\论文\201907TCDS\Code_and_Data\code\' testname '_mps.mat'])
load (['D:\OneDrive - mail.nankai.edu.cn\论文\201907TCDS\Code_and_Data\code\' testname '_mf.mat'])
%windows_number每个动作中的窗口数量
w_num=141;
%experiments_number实验次数
e_num=15;
%action_number动作数量
a_num=6;
%feature_demintion特征维数
f_dem=4;
time_domain_feature=[f_dem,w_num*e_num*a_num];
frequency_domain_feature=[f_dem,w_num*e_num*a_num];
time_frequency_domain_feature=[f_dem,w_num*e_num*a_num];
%正确率
accuracy_bp=zeros(c_t,7,3);
accuracy_ada=zeros(c_t,7,3);
%训练时间
train_time_bp=zeros(c_t,3);
train_time_ada=zeros(c_t,3);
%预测时间
test_time_bp=zeros(c_t,3);
test_time_ada=zeros(c_t,3);
%用于计算混淆矩阵
td_bp=zeros(c_t,635,6);
fd_bp=zeros(c_t,635,6);
tfd_bp=zeros(c_t,635,6);
td_ada=zeros(c_t,635,6);
fd_ada=zeros(c_t,635,6);
tfd_ada=zeros(c_t,635,6);
for exp_=1:15
    for act_=1:6
        %时域特征
        %rms第1通道
        time_domain_feature(1,(act_-1)*e_num*w_num+(exp_-1)*w_num+1:(act_-1)*e_num*w_num+exp_*w_num)=rms_feature(exp_,act_,:,1);
        %rms第2通道
        time_domain_feature(2,(act_-1)*e_num*w_num+(exp_-1)*w_num+1:(act_-1)*e_num*w_num+exp_*w_num)=rms_feature(exp_,act_,:,2);
        %zc第1通道
        time_domain_feature(3,(act_-1)*e_num*w_num+(exp_-1)*w_num+1:(act_-1)*e_num*w_num+exp_*w_num)=zc_feature(exp_,act_,:,1);
        %zc第2通道
        time_domain_feature(4,(act_-1)*e_num*w_num+(exp_-1)*w_num+1:(act_-1)*e_num*w_num+exp_*w_num)=zc_feature(exp_,act_,:,2);

        %频域特征
        %mps第1通道
        frequency_domain_feature(1,(act_-1)*e_num*w_num+(exp_-1)*w_num+1:(act_-1)*e_num*w_num+exp_*w_num)=mps_feature(exp_,act_,:,1);
        %mps第2通道
        frequency_domain_feature(2,(act_-1)*e_num*w_num+(exp_-1)*w_num+1:(act_-1)*e_num*w_num+exp_*w_num)=mps_feature(exp_,act_,:,2);
        %mf第1通道
        frequency_domain_feature(3,(act_-1)*e_num*w_num+(exp_-1)*w_num+1:(act_-1)*e_num*w_num+exp_*w_num)=mf_feature(exp_,act_,:,1);
        %mf第2通道
        frequency_domain_feature(4,(act_-1)*e_num*w_num+(exp_-1)*w_num+1:(act_-1)*e_num*w_num+exp_*w_num)=mf_feature(exp_,act_,:,2);
        
        %时域和频域特征
        %rms第1通道
        time_frequency_domain_feature(1,(act_-1)*e_num*w_num+(exp_-1)*w_num+1:(act_-1)*e_num*w_num+exp_*w_num)=rms_feature(exp_,act_,:,1);
        %rms第2通道
        time_frequency_domain_feature(2,(act_-1)*e_num*w_num+(exp_-1)*w_num+1:(act_-1)*e_num*w_num+exp_*w_num)=rms_feature(exp_,act_,:,2);
        %mps第1通道
        time_frequency_domain_feature(3,(act_-1)*e_num*w_num+(exp_-1)*w_num+1:(act_-1)*e_num*w_num+exp_*w_num)=mps_feature(exp_,act_,:,1);
        %mps第2通道
        time_frequency_domain_feature(4,(act_-1)*e_num*w_num+(exp_-1)*w_num+1:(act_-1)*e_num*w_num+exp_*w_num)=mps_feature(exp_,act_,:,2);
    end
end


%% 开始分类时域
for count_=1:c_t
    %% 生成分类标签
    label_=zeros(6,w_num*e_num*a_num);
    for i=1:6
        label_(i,(i-1)*w_num*e_num+1:(i)*w_num*e_num)=1;
    end
    %% 抽取训练集和测试集
    %将特征和标签合成一组
    feature_label=[time_domain_feature;label_];
    %训练比
    train_ratio=0.7;
    %用于训练的每个动作的数据数量
    train_number_action=fix(e_num*w_num*train_ratio);
    %用于测试的每个动作的数据数量
    test_number_action=e_num*w_num-train_number_action;
    %用于训练的数据总数量
    train_number=a_num*train_number_action;
    %用于测试的数据总数量
    test_number=test_number_action*a_num;
    %训练集
    train_set=zeros(a_num+f_dem,train_number);
    %测试集
    test_set=zeros(a_num+f_dem,test_number);
    for i=1:a_num%动作数量
        k=rand(e_num*w_num,1);
        [m,n]=sort(k);
        train_set(:,(i-1)*train_number_action+1:(i)*train_number_action)=feature_label(:,(i-1)*e_num*w_num+n(1:train_number_action));
        test_set(:,(i-1)*test_number_action+1:(i)*test_number_action)=feature_label(:,(i-1)*e_num*w_num+n(train_number_action+1:end));
    end
    %将训练集和测试集分为特征和标签
    train_feature=train_set(1:f_dem,:);
    train_label=train_set(f_dem+1:a_num+f_dem,:);
    test_feature=test_set(1:f_dem,:);
    test_label=test_set(f_dem+1:a_num+f_dem,:);
    train_classes=vec2ind(train_label);
    test_classes=vec2ind(test_label);
    %% 神经网络
    %bp训练时间
    tic
    net = patternnet(12);
    %修改至效果较好的训练函数
    net.trainFcn = 'trainlm';
    net.trainParam.showWindow=0; 
    %训练。matlab会自动对数据进行归一化，所以无需手动进行。
    net = train(net,train_feature,train_label);
    %仿真
    y=sim(net,test_feature);
    train_time_bp(count_,1)=toc;
    %bp预测时间
    tic
    %将输出归类
    predict_class_bp = vec2ind(y);
    td_bp(count_,:,:)=reshape(predict_class_bp,635,6);
    %计算准确率
    wrong_array_bp=(test_classes~=predict_class_bp);
    for i=1:6
        accuracy_bp(count_,i,1)=1-sum(wrong_array_bp((i-1)*test_number_action+1:i*test_number_action))/test_number_action;
    end
    accuracy_bp(count_,7,1)=mean(accuracy_bp(count_,1:6,1));
    test_time_bp(count_,1)=toc;
    % 显示进度
    clc
    disp (['BP时域循环次数：' num2str(count_)]);
    %% Adaboost
    %Adaboost训练时间
    tic
    ada=fitensemble(train_feature',train_classes,'AdaBoostM2',500,'Tree');
    train_time_ada(count_,1)=toc;
    tic
    predict_class_ada=predict(ada,test_feature');
    td_ada(count_,:,:)=reshape(predict_class_ada,635,6);
    %计算准确率
    wrong_array_ada=(test_classes~=predict_class_ada');
    for i=1:6
        accuracy_ada(count_,i,1)=1-sum(wrong_array_ada((i-1)*test_number_action+1:i*test_number_action))/test_number_action;
    end
    accuracy_ada(count_,7,1)=mean(accuracy_ada(count_,1:6,1));
    test_time_ada(count_,1)=toc;
    % 显示进度
    clc
    disp (['Adaboost时域循环次数：' num2str(count_)]);
end




%% 开始分类频域
for count_=1:c_t
    %% 生成分类标签
    label_=zeros(6,w_num*e_num*a_num);
    for i=1:6
        label_(i,(i-1)*w_num*e_num+1:(i)*w_num*e_num)=1;
    end
    %% 抽取训练集和测试集
    %将特征和标签合成一组
    feature_label=[frequency_domain_feature;label_];
    %训练比
    train_ratio=0.7;
    %用于训练的每个动作的数据数量
    train_number_action=fix(e_num*w_num*train_ratio);
    %用于测试的每个动作的数据数量
    test_number_action=e_num*w_num-train_number_action;
    %用于训练的数据总数量
    train_number=a_num*train_number_action;
    %用于测试的数据总数量
    test_number=test_number_action*a_num;
    %训练集
    train_set=zeros(a_num+f_dem,train_number);
    %测试集
    test_set=zeros(a_num+f_dem,test_number);
    for i=1:a_num%动作数量
        k=rand(e_num*w_num,1);
        [m,n]=sort(k);
        train_set(:,(i-1)*train_number_action+1:(i)*train_number_action)=feature_label(:,(i-1)*e_num*w_num+n(1:train_number_action));
        test_set(:,(i-1)*test_number_action+1:(i)*test_number_action)=feature_label(:,(i-1)*e_num*w_num+n(train_number_action+1:end));
    end
    %将训练集和测试集分为特征和标签
    train_feature=train_set(1:f_dem,:);
    train_label=train_set(f_dem+1:a_num+f_dem,:);
    test_feature=test_set(1:f_dem,:);
    test_label=test_set(f_dem+1:a_num+f_dem,:);
    train_classes=vec2ind(train_label);
    test_classes=vec2ind(test_label);
    %% 神经网络
    %bp训练时间
    tic
    net = patternnet(12);
    %修改至效果较好的训练函数
    net.trainFcn = 'trainlm';
    net.trainParam.showWindow=0; 
    %训练。matlab会自动对数据进行归一化，所以无需手动进行。
    net = train(net,train_feature,train_label);
    %仿真
    y=sim(net,test_feature);
    train_time_bp(count_,2)=toc;
    %bp预测时间
    tic
    %将输出归类
    predict_class_bp = vec2ind(y);
    fd_bp(count_,:,:)=reshape(predict_class_bp,635,6);
    %计算准确率
    wrong_array_bp=(test_classes~=predict_class_bp);
    for i=1:6
        accuracy_bp(count_,i,2)=1-sum(wrong_array_bp((i-1)*test_number_action+1:i*test_number_action))/test_number_action;
    end
    accuracy_bp(count_,7,2)=mean(accuracy_bp(count_,1:6,2));
    test_time_bp(count_,2)=toc;
    % 显示进度
    clc
    disp (['BP频域循环次数：' num2str(count_)]);
    %% Adaboost
    %Adaboost训练时间
    tic
    ada=fitensemble(train_feature',train_classes,'AdaBoostM2',500,'Tree');
    train_time_ada(count_,2)=toc;
    tic
    predict_class_ada=predict(ada,test_feature');
    fd_ada(count_,:,:)=reshape(predict_class_ada,635,6);
    %计算准确率
    wrong_array_ada=(test_classes~=predict_class_ada');
    for i=1:6
        accuracy_ada(count_,i,2)=1-sum(wrong_array_ada((i-1)*test_number_action+1:i*test_number_action))/test_number_action;
    end
    accuracy_ada(count_,7,2)=mean(accuracy_ada(count_,1:6,2));
    test_time_ada(count_,2)=toc;
    % 显示进度
    clc
    disp (['Adaboost频域循环次数：' num2str(count_)]);
end


%% 开始分类时频域
for count_=1:c_t
    %% 生成分类标签
    label_=zeros(6,w_num*e_num*a_num);
    for i=1:6
        label_(i,(i-1)*w_num*e_num+1:(i)*w_num*e_num)=1;
    end
    %% 抽取训练集和测试集
    %将特征和标签合成一组
    feature_label=[time_frequency_domain_feature;label_];
    %训练比
    train_ratio=0.7;
    %用于训练的每个动作的数据数量
    train_number_action=fix(e_num*w_num*train_ratio);
    %用于测试的每个动作的数据数量
    test_number_action=e_num*w_num-train_number_action;
    %用于训练的数据总数量
    train_number=a_num*train_number_action;
    %用于测试的数据总数量
    test_number=test_number_action*a_num;
    %训练集
    train_set=zeros(a_num+f_dem,train_number);
    %测试集
    test_set=zeros(a_num+f_dem,test_number);
    for i=1:a_num%动作数量
        k=rand(e_num*w_num,1);
        [m,n]=sort(k);
        train_set(:,(i-1)*train_number_action+1:(i)*train_number_action)=feature_label(:,(i-1)*e_num*w_num+n(1:train_number_action));
        test_set(:,(i-1)*test_number_action+1:(i)*test_number_action)=feature_label(:,(i-1)*e_num*w_num+n(train_number_action+1:end));
    end
    %将训练集和测试集分为特征和标签
    train_feature=train_set(1:f_dem,:);
    train_label=train_set(f_dem+1:a_num+f_dem,:);
    test_feature=test_set(1:f_dem,:);
    test_label=test_set(f_dem+1:a_num+f_dem,:);
    train_classes=vec2ind(train_label);
    test_classes=vec2ind(test_label);
    %% 神经网络
    %bp训练时间
    tic
    net = patternnet(12);
    %修改至效果较好的训练函数
    net.trainFcn = 'trainlm';
    net.trainParam.showWindow=0; 
    %训练。matlab会自动对数据进行归一化，所以无需手动进行。
    net = train(net,train_feature,train_label);
    %仿真
    y=sim(net,test_feature);
    train_time_bp(count_,3)=toc;
    %bp预测时间
    tic
    %将输出归类
    predict_class_bp = vec2ind(y);
    tfd_bp(count_,:,:)=reshape(predict_class_bp,635,6);
    %计算准确率
    wrong_array_bp=(test_classes~=predict_class_bp);
    for i=1:6
        accuracy_bp(count_,i,3)=1-sum(wrong_array_bp((i-1)*test_number_action+1:i*test_number_action))/test_number_action;
    end
    accuracy_bp(count_,7,3)=mean(accuracy_bp(count_,1:6,3));
    test_time_bp(count_,3)=toc;
    % 显示进度
    clc
    disp (['BP时频域循环次数：' num2str(count_)]);
    %% Adaboost
    %Adaboost训练时间
    tic
    ada=fitensemble(train_feature',train_classes,'AdaBoostM2',500,'Tree');
    train_time_ada(count_,3)=toc;
    tic
    predict_class_ada=predict(ada,test_feature');
    tfd_ada(count_,:,:)=reshape(predict_class_ada,635,6);
    %计算准确率
    wrong_array_ada=(test_classes~=predict_class_ada');
    for i=1:6
        accuracy_ada(count_,i,3)=1-sum(wrong_array_ada((i-1)*test_number_action+1:i*test_number_action))/test_number_action;
    end
    accuracy_ada(count_,7,3)=mean(accuracy_ada(count_,1:6,3));
    test_time_ada(count_,3)=toc;
    % 显示进度
    clc
    disp (['Adaboost时频域循环次数：' num2str(count_)]);
end



save([testname '_accuracy_bp.mat'],'accuracy_bp')
save([testname '_train_time_bp.mat'],'train_time_bp')
save([testname '_test_time_bp.mat'],'test_time_bp')
save([testname '_td_bp.mat'],'td_bp')
save([testname '_fd_bp.mat'],'fd_bp')
save([testname '_tfd_bp.mat'],'tfd_bp')

save([testname '_accuracy_ada.mat'],'accuracy_ada')
save([testname '_train_time_ada.mat'],'train_time_ada')
save([testname '_test_time_ada.mat'],'test_time_ada')
save([testname '_td_ada.mat'],'td_ada')
save([testname '_fd_ada.mat'],'fd_ada')
save([testname '_tfd_ada.mat'],'tfd_ada')
disp('all done!')