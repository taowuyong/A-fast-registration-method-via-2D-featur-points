%%A fast registration method via 2D featur points
str='E:\compile document\matlab\data\Indoor and outdoor dataset\Park32-Park14.txt';
fid = fopen(str,'r');
D = textscan(fid, '%f%f%f%f');
fclose(fid);
Ttrue=[D{1} D{2} D{3} D{4}];
pcloud1=pcread('E:\compile document\matlab\data\Indoor and outdoor dataset\Park32.ply');
PC1=pcloud1.Location;
pcloud2=pcread('E:\compile document\matlab\data\Indoor and outdoor dataset\Park14.ply');
PC2=pcloud2.Location;
tic;
% plot3(PC1(:,1),PC1(:,2),PC1(:,3),'.b','MarkerSize',1);
% hold on;
% plot3(PC2(:,1),PC2(:,2),PC2(:,3),'.g','MarkerSize',1);
% set(gca,'DataAspectRatio',[1 1 1]);
% axis off
% pr=0.022431135;          %apartment4  apartment6
% pr=0.025664071;          %boardroom0  boardroom3
% pr=0.16526765;          %castle1  castle2
% pr=0.12308180;          %City1  City2
pr=0.17055092;           %Park
[Mline1,Slinepoint1] = line2D( PC1,pr);
[Mline2,Slinepoint2] = line2D( PC2,pr);
% plot(Slinepoint1(:,1),Slinepoint1(:,2),'.b','MarkerSize',5);
% hold on;
% plot(Slinepoint2(:,1),Slinepoint2(:,2),'.g','MarkerSize',5);
% set(gca,'DataAspectRatio',[1 1 1]);
% axis off
[m1 n1]=size(Mline1);
MF2P1=[];
for i=1:m1-1
    for j=i+1:m1
       costheta=Mline1(i,1:2)*Mline1(j,1:2)';
       theta=acos(costheta)*180/pi;
       if theta>10 && theta<170
           a1=Mline1(i,1);
           b1=Mline1(i,2);
           c1=Mline1(i,3);
           a2=Mline1(j,1);
           b2=Mline1(j,2);
           c2=Mline1(j,3);
           x=(-c1*b2+b1*c2)/(a1*b2-b1*a2);
           y=(-a1*c2+c1*a2)/(a1*b2-b1*a2);
           MF2P1=[MF2P1;x y];
       end
    end
end
[m2 n2]=size(Mline2);
MF2P2=[];
for i=1:m2-1
    for j=i+1:m2
       costheta=Mline2(i,1:2)*Mline2(j,1:2)';
       theta=acos(costheta)*180/pi;
       if theta>10 && theta<170
           a1=Mline2(i,1);
           b1=Mline2(i,2);
           c1=Mline2(i,3);
           a2=Mline2(j,1);
           b2=Mline2(j,2);
           c2=Mline2(j,3);
           x=(-c1*b2+b1*c2)/(a1*b2-b1*a2);
           y=(-a1*c2+c1*a2)/(a1*b2-b1*a2);
           MF2P2=[MF2P2;x y];
       end
    end
end
% plot(Slinepoint1(:,1),Slinepoint1(:,2),'.b','MarkerSize',5);
% hold on;
% plot(MF2P1(:,1),MF2P1(:,2),'>r','MarkerSize',8);
% hold on;
% plot(Slinepoint2(:,1),Slinepoint2(:,2),'.g','MarkerSize',5);
% hold on;
% plot(MF2P2(:,1),MF2P2(:,2),'*k','MarkerSize',8);
% set(gca,'DataAspectRatio',[1 1 1]);
% axis off
sigma=3;    %%parameter, 1 for indoor point cloud and 3 for outdoor point cloud
[m3 n3]=size(MF2P1);
MFA1=[];
MFSL1=[];
for i=1:m3-2
    for j=i+1:m3-1
        L1=norm(MF2P1(i,:)-MF2P1(j,:));
        if L1>sigma     
            continue;
        end
        for k=j+1:m3
            L2=norm(MF2P1(j,:)-MF2P1(k,:));
            L3=norm(MF2P1(k,:)-MF2P1(i,:));
            if L2>sigma || L3>sigma      
                continue;
            end                             
            dL1=abs(L1-L2);
            dL2=abs(L2-L3);
            dL3=abs(L3-L2);
            if min([dL1 dL2 dL3])<3*pr
                continue;
            end
            L=[L1 L2 L3];
            idxv=[i j k];
            [Ls idxs]=sort(L);
            idxFA=[idxv(idxs(1)) idxv(idxs(2)) idxv(idxs(3))];
            MFA1=[MFA1;idxFA];
            MFSL1=[MFSL1;Ls];
        end
    end
end
[m4 n4]=size(MF2P2);
MFA2=[];
MFSL2=[];
for i=1:m4-2
    for j=i+1:m4-1
        L1=norm(MF2P2(i,:)-MF2P2(j,:));
        if L1>sigma  
            continue;
        end
        for k=j+1:m4
            L2=norm(MF2P2(j,:)-MF2P2(k,:));
            L3=norm(MF2P2(k,:)-MF2P2(i,:));
            if L2>sigma || L3>sigma    
                continue;
            end                              
            dL1=abs(L1-L2);
            dL2=abs(L2-L3);
            dL3=abs(L3-L2);
            if min([dL1 dL2 dL3])<3*pr
                continue;
            end
            L=[L1 L2 L3];
            idxv=[i j k];
            [Ls idxs]=sort(L);
            idxFA=[idxv(idxs(1)) idxv(idxs(2)) idxv(idxs(3))];
            MFA2=[MFA2;idxFA];
            MFSL2=[MFSL2;Ls];
        end
    end
end
m6=length(Slinepoint1);
m7=length(Slinepoint2);
[m5 n5]=size(MFA1);
[idx dist]=knnsearch(MFSL2,MFSL1,"K",1);  
overlap0=0;
for i=1:m5
    if dist(i)<0.2*pr      %%parameter
        FAS=MFA1(i,:);
        FAT=MFA2(idx(i),:);
        PbS=[MF2P1(FAS(1),:);MF2P1(FAS(2),:);MF2P1(FAS(3),:)]';
        PbT=[MF2P2(FAT(1),:);MF2P2(FAT(2),:);MF2P2(FAT(3),:)]';
        PmS=mean(PbS,2);
        PmT=mean(PbT,2);
        [U S V]=svd((PbS-PmS*[1 1 1])*(PbT-PmT*[1 1 1])');
        ra=V*U';
        ta=PmT-ra*PmS;
        tSlinepoint1=Slinepoint1*ra'+ones(m6,1)*ta';
        [idx1 dist1]=knnsearch(Slinepoint2,tSlinepoint1,"K",1);
        id=find(dist1<2*pr);
        overlap=length(id)/min(m6,m7);
        if overlap>overlap0
            r=ra;
            t=ta;
            overlap0=overlap;
        end
    end
end
tSlinepoint1=Slinepoint1*r'+ones(m6,1)*t';
% MF2P1t=MF2P1*r'+ones(m3,1)*t';
% plot(MF2P1t(:,1),MF2P1t(:,2),'>r','MarkerSize',8);
% hold on;
% plot(MF2P2(:,1),MF2P2(:,2),'ok','MarkerSize',8);
% hold on;
% plot(tSlinepoint1(:,1),tSlinepoint1(:,2),'.b','MarkerSize',5);
% hold on;
% plot(Slinepoint2(:,1),Slinepoint2(:,2),'.g','MarkerSize',5);
% hold on;
% set(gca,'DataAspectRatio',[1 1 1]);
% axis off
[idx2 dist2]=knnsearch(Slinepoint2,tSlinepoint1,'k',1);
idoverlap1=find(dist2<pr);
overlappoint=tSlinepoint1(idoverlap1,:);
Mtz=[];
for j=1:50
    Sidx=randperm(length(overlappoint),1);
    [idxx1 distt1]=rangesearch(PC1(:,1:2),Slinepoint1(idoverlap1(Sidx),:),2*pr);
    KNN1=PC1(idxx1{1},:);
    [idxx2 distt2]=rangesearch(PC2(:,1:2),overlappoint(Sidx,:),2*pr);
    KNN2=PC2(idxx2{1},:);
    zmin1=min(KNN1(:,3));
    zmin2=min(KNN2(:,3));
    tz=zmin2-zmin1;
    Mtz=[Mtz;tz 0];
end
MStz=[];
for i=1:1000
    seed=Mtz(1,:);
    Mtz(1,:)=[];
    for j=1:1000
        [idxtz disttz]=knnsearch(Mtz,seed,'k',1);
        [mi id]=min(disttz);
        if mi<0.2*pr         %参数
            seed=[seed;Mtz(idxtz(id),:)];
            Mtz(idxtz(id),:)=[];
        else
            break;
        end
    end
    [h l]=size(seed);
    Stz=[mean(seed(:,1)) h];
    MStz=[MStz;Stz];
    if isempty(Mtz)
        break
    end
end
[ma idStz]=max(MStz(:,2));
tz=MStz(idStz,1);
R=[r zeros(2,1);zeros(1,2) 1];
T=[t;tz];
TT=[R T;zeros(1,3) 1];
PC1t=PC1*R'+ones(length(PC1),1)*T';
time=toc;
plot3(PC1t(:,1),PC1t(:,2),PC1t(:,3),'.b','MarkerSize',1);
hold on;
plot3(PC2(:,1),PC2(:,2),PC2(:,3),'.g','MarkerSize',1);
set(gca,'DataAspectRatio',[1 1 1]);
axis off
Rtrue=Ttrue(1:3,1:3);
ttrue=Ttrue(1:3,4);
errorR=real(acos((trace(Rtrue*inv(R))-1)/2)*(180/pi))
errorH=norm(ttrue(1:2)-t)
errorV=norm(ttrue(3)-tz)
time