clc
clear
testname='yyk'

load (['D:\nutstore\TCDS_Special_Issue\Code_and_Data\features_rms\' testname '_rms.mat'])
for experiment_=1:15
    for action_=1:6
        x=reshape(rms_feature(experiment_,action_,:,1),[1 141]);
        y=reshape(rms_feature(experiment_,action_,:,2),[1 141]);
        mean_x(experiment_,action_)=mean(x);
        mean_y(experiment_,action_)=mean(y);
        r_rms(experiment_,action_)=mean(sqrt(x-mean(x)).^2+(y-mean(y)).^2);
    end
end
for experiment_=1:15
    for action_=1:6
        R_rms(experiment_,action_)=sum(sqrt((mean_x(experiment_,:)-mean_x(experiment_,action_)).^2+...
            (mean_y(experiment_,:)-mean_y(experiment_,action_)).^2))/5;
    end
end
sep_rms=mean(mean(r_rms))/mean(mean(R_rms))

load (['D:\nutstore\TCDS_Special_Issue\Code_and_Data\features_zc\' testname '_zc.mat'])
for experiment_=1:15
    for action_=1:6
        x=reshape(zc_feature(experiment_,action_,:,1),[1 141]);
        y=reshape(zc_feature(experiment_,action_,:,2),[1 141]);
        mean_x(experiment_,action_)=mean(x);
        mean_y(experiment_,action_)=mean(y);
        r_zc(experiment_,action_)=mean(sqrt(x-mean(x)).^2+(y-mean(y)).^2);
    end
end
for experiment_=1:15
    for action_=1:6
        R_zc(experiment_,action_)=sum(sqrt((mean_x(experiment_,:)-mean_x(experiment_,action_)).^2+...
            (mean_y(experiment_,:)-mean_y(experiment_,action_)).^2))/5;
    end
end
sep_zc=mean(mean(r_zc))/mean(mean(R_zc))

load (['D:\nutstore\TCDS_Special_Issue\Code_and_Data\features_wamp\' testname '_wamp.mat'])
for experiment_=1:15
    for action_=1:6
        x=reshape(wamp_feature(experiment_,action_,:,1),[1 141]);
        y=reshape(wamp_feature(experiment_,action_,:,2),[1 141]);
        mean_x(experiment_,action_)=mean(x);
        mean_y(experiment_,action_)=mean(y);
        r_wamp(experiment_,action_)=mean(sqrt(x-mean(x)).^2+(y-mean(y)).^2);
    end
end
for experiment_=1:15
    for action_=1:6
        R_wamp(experiment_,action_)=sum(sqrt((mean_x(experiment_,:)-mean_x(experiment_,action_)).^2+...
            (mean_y(experiment_,:)-mean_y(experiment_,action_)).^2))/5;
    end
end
sep_wamp=mean(mean(r_wamp))/mean(mean(R_wamp))

load (['D:\nutstore\TCDS_Special_Issue\Code_and_Data\features_mps\' testname '_mps.mat'])
for experiment_=1:15
    for action_=1:6
        x=reshape(mps_feature(experiment_,action_,:,1),[1 141]);
        y=reshape(mps_feature(experiment_,action_,:,2),[1 141]);
        mean_x(experiment_,action_)=mean(x);
        mean_y(experiment_,action_)=mean(y);
        r_mps(experiment_,action_)=mean(sqrt(x-mean(x)).^2+(y-mean(y)).^2);
    end
end
for experiment_=1:15
    for action_=1:6
        R_mps(experiment_,action_)=sum(sqrt((mean_x(experiment_,:)-mean_x(experiment_,action_)).^2+...
            (mean_y(experiment_,:)-mean_y(experiment_,action_)).^2))/5;
    end
end
sep_mps=mean(mean(r_mps))/mean(mean(R_mps))

load (['D:\nutstore\TCDS_Special_Issue\Code_and_Data\features_mf\' testname '_mf.mat'])
for experiment_=1:15
    for action_=1:6
        x=reshape(mf_feature(experiment_,action_,:,1),[1 141]);
        y=reshape(mf_feature(experiment_,action_,:,2),[1 141]);
        mean_x(experiment_,action_)=mean(x);
        mean_y(experiment_,action_)=mean(y);
        r_mf(experiment_,action_)=mean(sqrt(x-mean(x)).^2+(y-mean(y)).^2);
    end
end
for experiment_=1:15
    for action_=1:6
        R_mf(experiment_,action_)=sum(sqrt((mean_x(experiment_,:)-mean_x(experiment_,action_)).^2+...
            (mean_y(experiment_,:)-mean_y(experiment_,action_)).^2))/5;
    end
end
sep_mf=mean(mean(r_mf))/mean(mean(R_mf))

load (['D:\nutstore\TCDS_Special_Issue\Code_and_Data\features_mpf\' testname '_mpf.mat'])
for experiment_=1:15
    for action_=1:6
        x=reshape(mpf_feature(experiment_,action_,:,1),[1 141]);
        y=reshape(mpf_feature(experiment_,action_,:,2),[1 141]);
        mean_x(experiment_,action_)=mean(x);
        mean_y(experiment_,action_)=mean(y);
        r_mpf(experiment_,action_)=mean(sqrt(x-mean(x)).^2+(y-mean(y)).^2);
    end
end
for experiment_=1:15
    for action_=1:6
        R_mpf(experiment_,action_)=sum(sqrt((mean_x(experiment_,:)-mean_x(experiment_,action_)).^2+...
            (mean_y(experiment_,:)-mean_y(experiment_,action_)).^2))/5;
    end
end
sep_mpf=mean(mean(r_mpf))/mean(mean(R_mpf))