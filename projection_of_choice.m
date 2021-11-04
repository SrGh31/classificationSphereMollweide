%The hypersphere creation function is a slight variation of Roger Stafford's
%submission to MATLAB Central,-'Random points in an n-Dimensional
%Hypersphere', which can be found at
% https://nl.mathworks.com/matlabcentral/fileexchange/9443-random-points-in-an-n-dimensional-hypersphere?focused=5063757&tab=function
    
clear; clc;
addpath(genpath('~/PhD_work/Multi_verse/sc/angleGMLVQ/'));
addpath(genpath('~/PhD_work/Multi_verse/sc/fminlbfgs_version2c'));
load ('~/PhD_work/Multi_verse/sc/GCMS_expRatios.mat');
conditions = unique(GCMStab.CONDITION);
cls = zeros(size(GCMStab,1),1);
%load('Matrix_ALVQ_baseline_reged.mat'); % for the baseline version
load('CostFunct_Matrix.mat'); % for the cost weighted version
for j=1:length(conditions), cls(find(GCMStab.CONDITION==conditions(j))) = j;end
%finding indices of those samples which have labels either 1 or 7 or 8 or 9 only
 initIdx=find(cls==1 | cls==7 | cls==8 | cls==9);
    data_msng = allAnovaRatios(initIdx,:);
    label_select_cls = cls(initIdx);
    useIdx = find(sum(isnan(data_msng),2) < size(data_msng,2)*80/100);
    reqd_data=data_msng(useIdx,:); reqd_labels=label_select_cls(useIdx);
    idx=2;
    p_test=index_testing(:,idx); p_train=index_training(:,idx);
       p_test(p_test==0)=NaN; p_test=p_test(~isnan(p_test));
       p_train(p_train==0)=NaN; p_train=p_train(~isnan(p_train));
        train_label = reqd_labels(p_train); test_label = reqd_labels(p_test);
    prepro = struct('M',nanmean(reqd_data(p_train,:)),'SD',nanstd(reqd_data(p_train,:)));
    actAll = bsxfun(@rdivide,bsxfun(@minus,reqd_data,prepro.M),prepro.SD);
    actModel=MALVQ{23}; % for cw_ALVQ
    %actModel=MALVQ{4}; % for baseline
    actModel.omega=sqrt(actModel.A); 
    %actModel=Ensemble.expo;
    trainX=reqd_data(p_train,:);
   [estLabs,~,dist] = angleGMLVQ_classify(actAll,actModel);
   % [actModel,costs] = cw_angleGMLVQ_train(trainX,train_label,'beta',0.25,'dim',1,'PrototypesPerClass',1,'regularization',0,'costWeights',weights,'Display','off');
    x=actAll;
    x(isnan(x))=0;
    radius=1; nbr=size(x,1);
    %%
    x = bsxfun(@rdivide, (x*actModel.A'),sqrt(sum((x*actModel.A').^2)));
    s2 = sum(x.^2,2); 
    
    dimension=3;
    x= x.*repmat(radius*(ones(nbr,1).^(1/dimension))./sqrt(s2),1,dimension);
    prot=actModel.w;
    x_prot=bsxfun(@rdivide, (prot*actModel.A'),sqrt(sum((prot*actModel.A').^2)));
    s_prot = sum(x_prot.^2,2);
    x_prot= x_prot.*repmat(radius*(ones(4,1).^(1/dimension))./sqrt(s_prot),1,dimension);
t1=find(estLabs==1); t7=find(estLabs==7); t8=find(estLabs==8); t9=find(estLabs==9);
sz=[length(t1);length(t7);length(t8); length(t9)]; ind=NaN(sz(1),4);
ind(:,1)=t1; ind(1:sz(2),2)=t7; ind(1:sz(3),3)=t8; ind(1:sz(4),4)=t9; 
 prot_reduced=x_prot(:,1:3);
 % Dummy sphere creation
 samples=70000; dimension=3; label=[1,7,8,9];
     x2=rand(samples,dimension)*2-1; radius=1;
      sumsquared = @(X,dim) sum(X.^2, 2); prots=sumsquared(prot_reduced, 2);
      x2_sum=sumsquared(x2, 2);
       x2 = x2.*repmat(radius*(ones(samples,1).^(1/dimension))./sqrt(x2_sum),1,3);
       
for idk=1:samples
            dot_prod=x2(idk,:)*prot_reduced';
            ang_dist(idk,:)=bsxfun(@rdivide, dot_prod,sqrt(x2_sum(idk)*prots)');
            [dmin,pos]=max((ang_dist(idk,:)));
            estLabs2(idk)=label(pos);
end
%Setting perspectives
   A=90;B=-18;C=-26;
t21=find(estLabs2==1); t27=find(estLabs2==7); t28=find(estLabs2==8); t29=find(estLabs2==9);
sz=[length(t21);length(t27);length(t28); length(t29)]; ind2=NaN(sz(1),4);
ind2(:,1)=t21; ind2(1:sz(2),2)=t27; ind2(1:sz(3),3)=t28; ind2(1:sz(4),4)=t29; %col=['y','b','m','g'];
ind2(ind2==0)=NaN;
%Rotate viewpoint
[x2n,y2n,z2n]=perspective_change(x2(:,1),x2(:,2),x2(:,3),A,B,C);

%Option for projection
option=3;
%mapping out the dummy space
[azimuth_d,elevation_d,r_d] = cart2sph(x2n,y2n,z2n);
[X_d,Y_d]=map_projection(option,azimuth_d, elevation_d,r_d);
% mapping out real points and prototypes
[xn,yn,zn]=perspective_change(x(:,1),x(:,2),x(:,3),A,B,C);obj=['d','o','p','h'];
[azimuth_r,elevation_r,r_r] = cart2sph(xn,yn,zn);
[X_r,Y_r]=map_projection(option,azimuth_r, elevation_r,r_r);
[x_protn,y_protn,z_protn]=perspective_change(x_prot(:,1),x_prot(:,2),x_prot(:,3),A,B,C);
[azimuth_p,elevation_p,r_p] = cart2sph(x_protn,y_protn,z_protn);
[X_p,Y_p,str]=map_projection(option,azimuth_p, elevation_p,r_p);
%%
col={[0.3 0.7 0.7];[0.9 0.8 0.5];[0.7 0.6 0.9];[0.90 0.7 0.8]};

figure,
projections = 'perspectives_cw_ALVQ.gif';
pos1=[0 0.49 1 0.42];  ax(1)=subplot('Position',pos1);
col={[0.3 0.7 0.7];[0.9 0.8 0.5];[0.7 0.6 0.9];[0.90 0.7 0.8]};
for k = 1:length(sz) % loop over each non-empty color group
            indices=ind2(~isnan(ind2(:,k)),k);
            scatter(X_d(indices), Y_d(indices),5,col{k},...
                'Marker','s','MarkerFaceColor',col{k},'MarkerEdgeColor',col{k});
     hold on;
     %pause(2)
end
  % mapping out the real points
 
  count=1;col={[0.1055 0.6172 0.4648];[0.9 0.6 0.12];[0.4570 0.4375 0.6992];[0.9023 0.0602 0.5391]};
    
for k = 1:length(sz) % loop over each non-empty color group
            indices=ind(~isnan(ind(:,k)),k);
            scatter(X_r(indices), Y_r(indices),5,col{k},...
                'Marker','s','MarkerFaceColor',col{k},'MarkerEdgeColor',col{k});
     hold on;
     %pause(2)
end

for k = 1:length(sz) % loop over each non-empty color group
            %indices=ind2(~isnan(ind2(:,k)),k);
            scatter(X_p(k), Y_p(k),30,col{k},...
                'Marker',obj(k),'MarkerFaceColor',col{k},'MarkerEdgeColor','k');
     hold on;
     %pause(2)
end
axis square
 hold off;
 set(gca,'YTickLabel',[],'XTickLabel',[]);
 title([str,' perspective: ','X=',num2str(A),' Y=',num2str(B),' Z=', num2str(C)]);
 legend('dummy Healthy','dummy CYP21A2','dummy PORD','dummy SRD5A2','Healthy','Prot Healthy','CYP21A2','Prot CYP21A2','PORD','Prot PORD','SRD5A2','Prot SRD5A2','Location','BestOutside');
pos2=[0 0.01 0.95 0.42];  ax(2)=subplot('Position',pos2);
% dummy points on globe 
col={[0.3 0.7 0.7];[0.9 0.8 0.5];[0.7 0.6 0.9];[0.90 0.7 0.8]};
for k = 1:length(sz) % loop over each non-empty color group
            indices=ind2(~isnan(ind2(:,k)),k);
            scatter3(x2n(indices), y2n(indices),z2n(indices),5,col{k},...
                'Marker','s','MarkerFaceColor',col{k},'MarkerEdgeColor',col{k});
     hold on;
     %pause(2)
end
% real points and prototypes on globe 
col={[0.1055 0.6172 0.4648];[0.9 0.6 0.12];[0.4570 0.4375 0.6992];[0.9023 0.0602 0.5391]};

for k = 1:length(sz) % loop over each non-empty color group
            indices=ind(~isnan(ind(:,k)),k);
            scatter3(xn(indices), yn(indices),zn(indices),5,col{k},...
                'Marker','s','MarkerFaceColor',col{k},'MarkerEdgeColor',col{k});
     hold on;
     plot3(x_protn(k,1),y_protn(k,2),z_protn(k,3),'Marker',obj(k),'MarkerSize',10,'MarkerFaceColor','k','MarkerEdgeColor','y');
     %pause(2)
end
axis square
set(gca,'YTickLabel',[],'XTickLabel',[],'ZTickLabel',[]);
title(['Spherical perspective: ','X=',num2str(A),' Y=',num2str(B),' Z=', num2str(C)]);
set(get(gca,'title'),'Position',[0 -0.9 0.2 0.1]);
 legend('dummy Healthy','dummy CYP21A2','dummy PORD','dummy SRD5A2','Healthy','Prot Healthy','CYP21A2','Prot CYP21A2','PORD','Prot PORD','SRD5A2','Prot SRD5A2','Location','BestOutside');
%figure,

%linkprop(ax(1), 'CameraPosition');
pause(3)
%%
 for v=1:10:360         
            %view([v -15.2])
           % view([-10 v])
            view([65 v]) % for cw_ALVQ
            %view([v 190])  % for regularized baseline
            drawnow
            frame = getframe(1);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            %filename=sprintf(['PhD_work/Multi_verse/Animation/','filename_',num2str(count),'.png']);
            projections=sprintf(['PhD_work/Multi_verse/higherDimension/Visualization/','perspectives_cw_ALVQ.gif']);
            
            if v == 1
                  imwrite(imind,cm,projections,'gif','Loopcount',inf);
            else
                  imwrite(imind,cm,projections,'gif','WriteMode','append');
            end
           count=count+1;
            
end