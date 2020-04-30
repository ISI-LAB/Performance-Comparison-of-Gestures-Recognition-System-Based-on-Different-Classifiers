%���������Լ����BP�������Adaboost�ķ���������������ʱ��
clear;clc;close all;warning off;
%��������
testname='yyk';
%calculition_timeѭ������
c_t=1;
%��ȡ�����ǰ��շ���ȵĴ�С���������
load (['D:\OneDrive - mail.nankai.edu.cn\����\201907TCDS\Code_and_Data\code\' testname '_rms.mat'])
load (['D:\OneDrive - mail.nankai.edu.cn\����\201907TCDS\Code_and_Data\code\' testname '_zc.mat'])
load (['D:\OneDrive - mail.nankai.edu.cn\����\201907TCDS\Code_and_Data\code\' testname '_mps.mat'])
load (['D:\OneDrive - mail.nankai.edu.cn\����\201907TCDS\Code_and_Data\code\' testname '_mf.mat'])
%windows_numberÿ�������еĴ�������
w_num=141;
%experiments_numberʵ�����
e_num=15;
%action_number��������
a_num=6;
%feature_demintion����ά��
f_dem=4;
time_domain_feature=[f_dem,w_num*e_num*a_num];
frequency_domain_feature=[f_dem,w_num*e_num*a_num];
time_frequency_domain_feature=[f_dem,w_num*e_num*a_num];
%��ȷ��
accuracy_bp=zeros(c_t,7,3);
accuracy_ada=zeros(c_t,7,3);
%ѵ��ʱ��
train_time_bp=zeros(c_t,3);
train_time_ada=zeros(c_t,3);
%Ԥ��ʱ��
test_time_bp=zeros(c_t,3);
test_time_ada=zeros(c_t,3);
%���ڼ����������
td_bp=zeros(c_t,635,6);
fd_bp=zeros(c_t,635,6);
tfd_bp=zeros(c_t,635,6);
td_ada=zeros(c_t,635,6);
fd_ada=zeros(c_t,635,6);
tfd_ada=zeros(c_t,635,6);
for exp_=1:15
    for act_=1:6
        %ʱ������
        %rms��1ͨ��
        time_domain_feature(1,(act_-1)*e_num*w_num+(exp_-1)*w_num+1:(act_-1)*e_num*w_num+exp_*w_num)=rms_feature(exp_,act_,:,1);
        %rms��2ͨ��
        time_domain_feature(2,(act_-1)*e_num*w_num+(exp_-1)*w_num+1:(act_-1)*e_num*w_num+exp_*w_num)=rms_feature(exp_,act_,:,2);
        %zc��1ͨ��
        time_domain_feature(3,(act_-1)*e_num*w_num+(exp_-1)*w_num+1:(act_-1)*e_num*w_num+exp_*w_num)=zc_feature(exp_,act_,:,1);
        %zc��2ͨ��
        time_domain_feature(4,(act_-1)*e_num*w_num+(exp_-1)*w_num+1:(act_-1)*e_num*w_num+exp_*w_num)=zc_feature(exp_,act_,:,2);

        %Ƶ������
        %mps��1ͨ��
        frequency_domain_feature(1,(act_-1)*e_num*w_num+(exp_-1)*w_num+1:(act_-1)*e_num*w_num+exp_*w_num)=mps_feature(exp_,act_,:,1);
        %mps��2ͨ��
        frequency_domain_feature(2,(act_-1)*e_num*w_num+(exp_-1)*w_num+1:(act_-1)*e_num*w_num+exp_*w_num)=mps_feature(exp_,act_,:,2);
        %mf��1ͨ��
        frequency_domain_feature(3,(act_-1)*e_num*w_num+(exp_-1)*w_num+1:(act_-1)*e_num*w_num+exp_*w_num)=mf_feature(exp_,act_,:,1);
        %mf��2ͨ��
        frequency_domain_feature(4,(act_-1)*e_num*w_num+(exp_-1)*w_num+1:(act_-1)*e_num*w_num+exp_*w_num)=mf_feature(exp_,act_,:,2);
        
        %ʱ���Ƶ������
        %rms��1ͨ��
        time_frequency_domain_feature(1,(act_-1)*e_num*w_num+(exp_-1)*w_num+1:(act_-1)*e_num*w_num+exp_*w_num)=rms_feature(exp_,act_,:,1);
        %rms��2ͨ��
        time_frequency_domain_feature(2,(act_-1)*e_num*w_num+(exp_-1)*w_num+1:(act_-1)*e_num*w_num+exp_*w_num)=rms_feature(exp_,act_,:,2);
        %mps��1ͨ��
        time_frequency_domain_feature(3,(act_-1)*e_num*w_num+(exp_-1)*w_num+1:(act_-1)*e_num*w_num+exp_*w_num)=mps_feature(exp_,act_,:,1);
        %mps��2ͨ��
        time_frequency_domain_feature(4,(act_-1)*e_num*w_num+(exp_-1)*w_num+1:(act_-1)*e_num*w_num+exp_*w_num)=mps_feature(exp_,act_,:,2);
    end
end


%% ��ʼ����ʱ��
for count_=1:c_t
    %% ���ɷ����ǩ
    label_=zeros(6,w_num*e_num*a_num);
    for i=1:6
        label_(i,(i-1)*w_num*e_num+1:(i)*w_num*e_num)=1;
    end
    %% ��ȡѵ�����Ͳ��Լ�
    %�������ͱ�ǩ�ϳ�һ��
    feature_label=[time_domain_feature;label_];
    %ѵ����
    train_ratio=0.7;
    %����ѵ����ÿ����������������
    train_number_action=fix(e_num*w_num*train_ratio);
    %���ڲ��Ե�ÿ����������������
    test_number_action=e_num*w_num-train_number_action;
    %����ѵ��������������
    train_number=a_num*train_number_action;
    %���ڲ��Ե�����������
    test_number=test_number_action*a_num;
    %ѵ����
    train_set=zeros(a_num+f_dem,train_number);
    %���Լ�
    test_set=zeros(a_num+f_dem,test_number);
    for i=1:a_num%��������
        k=rand(e_num*w_num,1);
        [m,n]=sort(k);
        train_set(:,(i-1)*train_number_action+1:(i)*train_number_action)=feature_label(:,(i-1)*e_num*w_num+n(1:train_number_action));
        test_set(:,(i-1)*test_number_action+1:(i)*test_number_action)=feature_label(:,(i-1)*e_num*w_num+n(train_number_action+1:end));
    end
    %��ѵ�����Ͳ��Լ���Ϊ�����ͱ�ǩ
    train_feature=train_set(1:f_dem,:);
    train_label=train_set(f_dem+1:a_num+f_dem,:);
    test_feature=test_set(1:f_dem,:);
    test_label=test_set(f_dem+1:a_num+f_dem,:);
    train_classes=vec2ind(train_label);
    test_classes=vec2ind(test_label);
    %% ������
    %bpѵ��ʱ��
    tic
    net = patternnet(12);
    %�޸���Ч���Ϻõ�ѵ������
    net.trainFcn = 'trainlm';
    net.trainParam.showWindow=0; 
    %ѵ����matlab���Զ������ݽ��й�һ�������������ֶ����С�
    net = train(net,train_feature,train_label);
    %����
    y=sim(net,test_feature);
    train_time_bp(count_,1)=toc;
    %bpԤ��ʱ��
    tic
    %���������
    predict_class_bp = vec2ind(y);
    td_bp(count_,:,:)=reshape(predict_class_bp,635,6);
    %����׼ȷ��
    wrong_array_bp=(test_classes~=predict_class_bp);
    for i=1:6
        accuracy_bp(count_,i,1)=1-sum(wrong_array_bp((i-1)*test_number_action+1:i*test_number_action))/test_number_action;
    end
    accuracy_bp(count_,7,1)=mean(accuracy_bp(count_,1:6,1));
    test_time_bp(count_,1)=toc;
    % ��ʾ����
    clc
    disp (['BPʱ��ѭ��������' num2str(count_)]);
    %% Adaboost
    %Adaboostѵ��ʱ��
    tic
    ada=fitensemble(train_feature',train_classes,'AdaBoostM2',500,'Tree');
    train_time_ada(count_,1)=toc;
    tic
    predict_class_ada=predict(ada,test_feature');
    td_ada(count_,:,:)=reshape(predict_class_ada,635,6);
    %����׼ȷ��
    wrong_array_ada=(test_classes~=predict_class_ada');
    for i=1:6
        accuracy_ada(count_,i,1)=1-sum(wrong_array_ada((i-1)*test_number_action+1:i*test_number_action))/test_number_action;
    end
    accuracy_ada(count_,7,1)=mean(accuracy_ada(count_,1:6,1));
    test_time_ada(count_,1)=toc;
    % ��ʾ����
    clc
    disp (['Adaboostʱ��ѭ��������' num2str(count_)]);
end




%% ��ʼ����Ƶ��
for count_=1:c_t
    %% ���ɷ����ǩ
    label_=zeros(6,w_num*e_num*a_num);
    for i=1:6
        label_(i,(i-1)*w_num*e_num+1:(i)*w_num*e_num)=1;
    end
    %% ��ȡѵ�����Ͳ��Լ�
    %�������ͱ�ǩ�ϳ�һ��
    feature_label=[frequency_domain_feature;label_];
    %ѵ����
    train_ratio=0.7;
    %����ѵ����ÿ����������������
    train_number_action=fix(e_num*w_num*train_ratio);
    %���ڲ��Ե�ÿ����������������
    test_number_action=e_num*w_num-train_number_action;
    %����ѵ��������������
    train_number=a_num*train_number_action;
    %���ڲ��Ե�����������
    test_number=test_number_action*a_num;
    %ѵ����
    train_set=zeros(a_num+f_dem,train_number);
    %���Լ�
    test_set=zeros(a_num+f_dem,test_number);
    for i=1:a_num%��������
        k=rand(e_num*w_num,1);
        [m,n]=sort(k);
        train_set(:,(i-1)*train_number_action+1:(i)*train_number_action)=feature_label(:,(i-1)*e_num*w_num+n(1:train_number_action));
        test_set(:,(i-1)*test_number_action+1:(i)*test_number_action)=feature_label(:,(i-1)*e_num*w_num+n(train_number_action+1:end));
    end
    %��ѵ�����Ͳ��Լ���Ϊ�����ͱ�ǩ
    train_feature=train_set(1:f_dem,:);
    train_label=train_set(f_dem+1:a_num+f_dem,:);
    test_feature=test_set(1:f_dem,:);
    test_label=test_set(f_dem+1:a_num+f_dem,:);
    train_classes=vec2ind(train_label);
    test_classes=vec2ind(test_label);
    %% ������
    %bpѵ��ʱ��
    tic
    net = patternnet(12);
    %�޸���Ч���Ϻõ�ѵ������
    net.trainFcn = 'trainlm';
    net.trainParam.showWindow=0; 
    %ѵ����matlab���Զ������ݽ��й�һ�������������ֶ����С�
    net = train(net,train_feature,train_label);
    %����
    y=sim(net,test_feature);
    train_time_bp(count_,2)=toc;
    %bpԤ��ʱ��
    tic
    %���������
    predict_class_bp = vec2ind(y);
    fd_bp(count_,:,:)=reshape(predict_class_bp,635,6);
    %����׼ȷ��
    wrong_array_bp=(test_classes~=predict_class_bp);
    for i=1:6
        accuracy_bp(count_,i,2)=1-sum(wrong_array_bp((i-1)*test_number_action+1:i*test_number_action))/test_number_action;
    end
    accuracy_bp(count_,7,2)=mean(accuracy_bp(count_,1:6,2));
    test_time_bp(count_,2)=toc;
    % ��ʾ����
    clc
    disp (['BPƵ��ѭ��������' num2str(count_)]);
    %% Adaboost
    %Adaboostѵ��ʱ��
    tic
    ada=fitensemble(train_feature',train_classes,'AdaBoostM2',500,'Tree');
    train_time_ada(count_,2)=toc;
    tic
    predict_class_ada=predict(ada,test_feature');
    fd_ada(count_,:,:)=reshape(predict_class_ada,635,6);
    %����׼ȷ��
    wrong_array_ada=(test_classes~=predict_class_ada');
    for i=1:6
        accuracy_ada(count_,i,2)=1-sum(wrong_array_ada((i-1)*test_number_action+1:i*test_number_action))/test_number_action;
    end
    accuracy_ada(count_,7,2)=mean(accuracy_ada(count_,1:6,2));
    test_time_ada(count_,2)=toc;
    % ��ʾ����
    clc
    disp (['AdaboostƵ��ѭ��������' num2str(count_)]);
end


%% ��ʼ����ʱƵ��
for count_=1:c_t
    %% ���ɷ����ǩ
    label_=zeros(6,w_num*e_num*a_num);
    for i=1:6
        label_(i,(i-1)*w_num*e_num+1:(i)*w_num*e_num)=1;
    end
    %% ��ȡѵ�����Ͳ��Լ�
    %�������ͱ�ǩ�ϳ�һ��
    feature_label=[time_frequency_domain_feature;label_];
    %ѵ����
    train_ratio=0.7;
    %����ѵ����ÿ����������������
    train_number_action=fix(e_num*w_num*train_ratio);
    %���ڲ��Ե�ÿ����������������
    test_number_action=e_num*w_num-train_number_action;
    %����ѵ��������������
    train_number=a_num*train_number_action;
    %���ڲ��Ե�����������
    test_number=test_number_action*a_num;
    %ѵ����
    train_set=zeros(a_num+f_dem,train_number);
    %���Լ�
    test_set=zeros(a_num+f_dem,test_number);
    for i=1:a_num%��������
        k=rand(e_num*w_num,1);
        [m,n]=sort(k);
        train_set(:,(i-1)*train_number_action+1:(i)*train_number_action)=feature_label(:,(i-1)*e_num*w_num+n(1:train_number_action));
        test_set(:,(i-1)*test_number_action+1:(i)*test_number_action)=feature_label(:,(i-1)*e_num*w_num+n(train_number_action+1:end));
    end
    %��ѵ�����Ͳ��Լ���Ϊ�����ͱ�ǩ
    train_feature=train_set(1:f_dem,:);
    train_label=train_set(f_dem+1:a_num+f_dem,:);
    test_feature=test_set(1:f_dem,:);
    test_label=test_set(f_dem+1:a_num+f_dem,:);
    train_classes=vec2ind(train_label);
    test_classes=vec2ind(test_label);
    %% ������
    %bpѵ��ʱ��
    tic
    net = patternnet(12);
    %�޸���Ч���Ϻõ�ѵ������
    net.trainFcn = 'trainlm';
    net.trainParam.showWindow=0; 
    %ѵ����matlab���Զ������ݽ��й�һ�������������ֶ����С�
    net = train(net,train_feature,train_label);
    %����
    y=sim(net,test_feature);
    train_time_bp(count_,3)=toc;
    %bpԤ��ʱ��
    tic
    %���������
    predict_class_bp = vec2ind(y);
    tfd_bp(count_,:,:)=reshape(predict_class_bp,635,6);
    %����׼ȷ��
    wrong_array_bp=(test_classes~=predict_class_bp);
    for i=1:6
        accuracy_bp(count_,i,3)=1-sum(wrong_array_bp((i-1)*test_number_action+1:i*test_number_action))/test_number_action;
    end
    accuracy_bp(count_,7,3)=mean(accuracy_bp(count_,1:6,3));
    test_time_bp(count_,3)=toc;
    % ��ʾ����
    clc
    disp (['BPʱƵ��ѭ��������' num2str(count_)]);
    %% Adaboost
    %Adaboostѵ��ʱ��
    tic
    ada=fitensemble(train_feature',train_classes,'AdaBoostM2',500,'Tree');
    train_time_ada(count_,3)=toc;
    tic
    predict_class_ada=predict(ada,test_feature');
    tfd_ada(count_,:,:)=reshape(predict_class_ada,635,6);
    %����׼ȷ��
    wrong_array_ada=(test_classes~=predict_class_ada');
    for i=1:6
        accuracy_ada(count_,i,3)=1-sum(wrong_array_ada((i-1)*test_number_action+1:i*test_number_action))/test_number_action;
    end
    accuracy_ada(count_,7,3)=mean(accuracy_ada(count_,1:6,3));
    test_time_ada(count_,3)=toc;
    % ��ʾ����
    clc
    disp (['AdaboostʱƵ��ѭ��������' num2str(count_)]);
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