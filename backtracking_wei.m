

load 2017.LONGITUDE;               %不用变
load 2017.LATITUDE;                %不用变
load Depth;                        %不用变
load 2019.U_VELOCITY_3m.mat;
load 2019.V_VELOCITY_3m.mat;
load 2019.HEIGHT_3m.mat;           
load 2019.TEMPERATURE_bot.mat;     
load 2019_Genetic_info_3m;
load Grid_reshape_longitude;       %不用变
load Grid_reshape_latitude;        %不用变

%%实际水深
ADEPTH=Depth-HEIGHT;


%%网格                        
lon1d=LONGITUDE_reshape;      %reshape过的网格     366×156                             
lat1d=LATITUDE_reshape;


%% Get land mask for boundary conditions
%1 land, 0 water
%grid t u v have the same land mask
imask=zeros(length(LONGITUDE(:,1)),length(LONGITUDE(1,:)));                       % imask (366×802)
imask(isnan(LONGITUDE))=1; %land 1                                  % replace the NAN with 1 in LONGITUDE, and other elements with 0 in LONGITUDE,从而组成一个新的矩阵即imask。就是说imask中1为陆地，0为水
imask=imask(:,15:170);                                                               % imask (366×156)

imask1d1=imask(:);                                                       %按列，第一列第二列第三列。。。输入到in.imask1d1中
imask1d1(1)=1; %set imask1d(1)=1 as land, NaN will lead iw=1, so idxP=find(imask1dtmp(iw)==0) is empty


%% 取不同的case
aa=1;
bb=369;
casen=length(Genetic_info(aa:bb,3));     % case数量
jcount=0;
         
      pnumber=100;  % 初始投放粒子数
      kx=0.1;    %horizontal diffusivity
      ky=0.1;
 %对应178行

all_templon=zeros(casen,28*24*6,pnumber);     % itotalcount: itd*nsteps, itd=24*age, max(age)=42
all_templat=zeros(casen,28*24*6,pnumber);
% all_temph=zeros(casen,43*24*6,pnumber);
all_tempt=zeros(casen,28*24*6,pnumber);
all_tempah=zeros(casen,28*24*6,pnumber);
all_tempU=zeros(casen,28*24*6,pnumber);
all_tempV=zeros(casen,28*24*6,pnumber);

hatch_templon=zeros(casen,1,pnumber);                                     %每件case的最终的位置即hatch location
hatch_templat=zeros(casen,1,pnumber);
%hatch_temph=zeros(casen,1,pnumber);
hatch_tempt=zeros(casen,1,pnumber);
hatch_tempah=zeros(casen,1,pnumber);
hatch_tempU=zeros(casen,1,pnumber);
hatch_tempV=zeros(casen,1,pnumber);

% all_actlon=zeros(casen,1,pnumber);
% all_actlat=zeros(casen,1,pnumber);

age=Genetic_info(aa:bb,8);            %幼鱼的年纪
si=Genetic_info(aa:bb,9);             %模拟的初始时刻
initial_lon=Genetic_info(aa:bb,3);     %模拟的初始位置，即幼鱼的捕捉位置
initial_lat=Genetic_info(aa:bb,2);

for ii=1:1:casen
    jcount=jcount+1;
    dis_drift=0;
disp('now: ii')    % 显示实时的ii数
[ii]

u1=U_VELOCITY(:,15:170,1:si(jcount));                                         %U_VELOCITY (366×156×348...)
v1=V_VELOCITY(:,15:170,1:si(jcount));                                         
%h1=HEIGHT(:,15:170,1:si(jcount));
t1=TEMPERATURE(:,15:170,1:si(jcount));
ah1=ADEPTH(:,15:170,1:si(jcount));

u1=permute(u1, [3 1 2]);                                              %in.u1从366×802×348 变为348×366×802， 366×802为网格数，366为时间
v1=permute(v1, [3 1 2]);
%h1=permute(h1, [3 1 2]);
t1=permute(t1, [3 1 2]);
ah1=permute(ah1,[3 1 2]);

u1=flip(u1,1);                              % wei，将流速对时间转置
v1=flip(v1,1);
%h1=flip(h1,1);
t1=flip(t1,1);
ah1=flip(ah1,1);


      Lon1d=zeros(pnumber,1); 
      Lat1d=zeros(pnumber,1);

      Lon1d(1:pnumber,1)=initial_lon(jcount);  
      Lat1d(1:pnumber,1)=initial_lat(jcount);


    Lon = Lon1d;                                                                   % in.Lon1d和in.Lat1d为初始时刻投入的颗粒的位置,setting from the file "initParticle.m"
    Lat = Lat1d;                                                                   % Lon和Lat的值每一次循环之后会变化，代表瞬态位置，而in.Lon1d和in.Lat1d在本次运行中的值并不会变化，他们只会在再次运行之后变化，这是由 "initParticle.m"决定的。
                                                                                      % 因此，St_templon(1,:)=in.Lon1d, St_templon(in.nsteps*in.itd+1,:)=Lon,就是最后一行等于Lon，Lat同理。

      
 
      itd=24*age(jcount);     %模拟的时长
      nsteps = 6;   % Number of steps within a pentad
      dt=(1/nsteps)*60*60;        % 10 mins, time steps; 1hour: time interval of flow field output from AEM3D
      dthour=dt/3600;
      totalicount=1+nsteps*(itd);

 



	icount=1;                                                                             %这里icount应该为0吧？NO，应该是1！
%     St_driftx=zeros(totalicount,pnumber);
%     St_drifty=zeros(totalicount,pnumber);
%     St_rwalkx=zeros(totalicount,pnumber);
%     St_rwalky=zeros(totalicount,pnumber);

	St_templat=zeros(totalicount,pnumber);    
	St_templon=zeros(totalicount,pnumber);
	St_tempT=zeros(totalicount,pnumber);
	St_tempS=zeros(totalicount,pnumber);
	St_tempU=zeros(totalicount,pnumber);
	St_tempV=zeros(totalicount,pnumber);                                 %应该是三维的
%     St_tempH=zeros(totalicount,pnumber);
    St_tempAH=zeros(totalicount,pnumber);
 
    St_giw=zeros(totalicount,pnumber);


	St_templat(1,:)=Lat1d;
	St_templon(1,:)=Lon1d;
%totalicount:total time intervals，总的时间间隔，=2001（总模拟时长为250day，即6000h，dt=3h，故有6000/3+1=2001次循环）；
%pnumber：比如同时release500个颗粒，即为500，=5；St_templat即为初始时刻5个颗粒的纬度位置    
    

it1=1;it2=itd;

for it = it1:it2                                                             %it2-it1+1=in.itd+1=49+1=50                

   % pericount=100*(icount/totalicount);
   % disp('Start:it1 End:it2 now:it % ii')
    %[it1 it2 it pericount]

    for k = 1:nsteps                          
totalx=0;
totaly=0;

        icount=icount+1;
                varalpha = k / nsteps;                    % nsteps=40; varalpha<=1                                  %nsteps是啥？k是啥？ varalpha是啥？


        u2 = squeeze((1-varalpha) * u1(it,:,:) + varalpha * u1(it+1,:,:));                          % u1: 169×87×193, it: 37:86; u2: 87×193
        v2 = squeeze((1-varalpha) * v1(it,:,:) + varalpha * v1(it+1,:,:));                          % u1 or 流场的时间步长是3h*40=120h=5days
      %  h2 = squeeze((1-varalpha) * h1(it,:,:) + varalpha * h1(it+1,:,:));
        t2 = squeeze((1-varalpha) * t1(it,:,:) + varalpha * t1(it+1,:,:));
        ah2= squeeze((1-varalpha) * ah1(it,:,:) + varalpha * ah1(it+1,:,:));

%for Boundary conditon
	imask1dtmp=imask1d1;



%to speed up the program, using a subfield
%
subuffer=0.5;                                                                                                                         
minlon=min(Lon)-subuffer;
minlat=min(Lat)-subuffer;
maxlon=max(Lon)+subuffer;
maxlat=max(Lat)+subuffer;

isub=find(minlon<lon1d&lon1d<maxlon &  minlat<lat1d&lat1d<maxlat);              % isub是索引。为什么在（in.lon,in.lat）这个范围内找？因为如果用原来的网格（in.lon1d,in.lat1d），inlonsub/inlatsub会出现NaN值，用（in.lon,in.lat）就不会，且不会影响
                                                                                    % U/V的取值。因为(Lon, Lat)一直在原来的网格内，在水域内，那么离他最近的点也一定永远在水域内，U/V就永远不会是NaN值。但是若是particle运动到水域的边缘呢？离他最近的点一定还在水域内吗？
inlonsub=lon1d(isub);                                                     % 坐标值
inlatsub=lat1d(isub);
imask1dtmpsub=imask1dtmp(isub);



                                                                                                    % 原理是什么？为什么用距其距离最近的点，为什么不直接用该点的流速？如何保证最近的点位流速不是NAN？
if icount==2
iw=myfindgrid(Lon,Lat,inlonsub,inlatsub);% where the particle located at that time                     % 找出距离（五个）初始点位（Lon,Lat）距离最短的点!即距五个点最近的点，共五个
iwold1=iw;                                                                                             %iw即为最近的点在isub中的索引位置
giw=isub(iw);                                                                                          %isub中iw位置处的值为giw，该值代入in.lon×in.lat中即可得到最近点的坐标值，或者代入流场中得到最近点的流速值。
end

	%assign model grid velocity to each particle
        U=u2(giw);                                                                   %U=NAN, V=NAN, T=NAN，咋整？u2:87×193
        V=v2(giw);
      %  H=h2(giw);
        T=t2(giw);
        AH=ah2(giw);

        driftx=(V*dt)./(111000*cosd(Lat));                                             % V正方向是向东，即向右                % cosd(Lat)求Lat的cosine值  1经度的距离约等于111.11cos(纬度θ)km。driftx与drifty都是球坐标，即经纬度，这样才能与Lon,Lat相加！
        drifty=-(U*dt)./(111000);                                % U正方向是向南，即向下                %  1纬度的距离约等于111.11km;
        totalx=totalx-driftx;
        totaly=totaly-drifty;  
                                                                   
	rwalkx=randn(size(Lon))*sqrt(2*kx*dt)./(111000*cosd(Lat));
	%rwalky=(2.*rand(size(Lat))-1)*sqrt(2*ky*dt)/111000;
    rwalky=randn(size(Lat))*sqrt(2*ky*dt)./111000;
        totalx=totalx+rwalkx;                                        
        totaly=totaly+rwalky;                                        

% for ii=1:pnumber
% if isnan(V(ii))
%     V(ii)=0;U(ii)=0;
% end
% end


% A=V*dt;
% B=U*dt;
% drift=sqrt(A.^2+B.^2);
% dis_drift=drift+dis_drift;       %求粒子运动的轨迹的总距离



        iw=myfindgrid(Lon+totalx,Lat+totaly,inlonsub,inlatsub);                               
      
        idxP=find(imask1dtmpsub(iw)==0); %particle in water after the totalx/y added                           %1 land, 0 water
        %only move partilce to water point
       
        Lon(idxP)=Lon(idxP)+totalx(idxP);
        Lat(idxP)=Lat(idxP)+totaly(idxP);                %一旦颗粒物进入陆地，颗粒物保持不动等待下一个循环继续运动

%Boundary2 一旦释放的100颗粒全都进入陆地，结束该次循环；如果只是部分颗粒进入陆地，让其经纬度变为nan值，其余颗粒继续循环。       
        idxP2=find(imask1dtmpsub(iw)==1);
         Lon(idxP2)=NaN;
         Lat(idxP2)=NaN;
    if (imask1dtmpsub(iw)==1)  
        break;
    end
        
       

        iwold2=iw;
 if icount==2
        iwold=iwold1;
 end
        iwold(idxP)=iwold2(idxP);
    
        %convert iwold from sub to global
        giw(idxP)=isub(iwold(idxP));






%St_idxP(icount,:)=imask1dtmpsub(iw);


St_giw(icount,:)=giw(:);
St_templat(icount,:)=Lat(:);
St_templon(icount,:)=Lon(:);                     % St_templon(2,1)-St_templon(1,1)=St_totalx(2,1)
%St_tempT(icount,:)=T;
%St_tempS(icount,:)=S;
St_tempU(icount,:)=U;
St_tempV(icount,:)=V;
%St_tempH(icount,:)=H;
St_tempT(icount,:)=T;
St_tempAH(icount,:)=AH;

% St_driftx(icount,:)=driftx(:);
% St_drifty(icount,:)=drifty(:);
% St_rwalkx(icount,:)=rwalkx(:);
% St_rwalky(icount,:)=rwalky(:);




    end %nsteps

     if (imask1dtmpsub(iw)==1)
        break;
    end


end %it1

% iicount=0;
% for ii=1:pnumber
%     iicount=1+iicount;
%     ax=St_driftx(:,iicount);  ay=St_drifty(:,iicount);
%     arx=St_rwalkx(:,iicount); ary=St_rwalky(:,iicount);
%     ax(isnan(ax))=[];   ay(isnan(ay))=[];
%     arx(isnan(arx))=[]; ary(isnan(ary))=[];
%     [m n]=size(ax);
%     bx=ax(m);   by=ay(m);
%     brx=arx(m); bry=ary(m);
%     cx(1,iicount)=bx;   cy(1,iicount)=by;
%     crx(1,iicount)=brx; cry(1,iicount)=bry; 
% end
% act_lon=St_templon(icount,:)+cx-crx;
% act_lat=St_templat(icount,:)+cy-cry;
% 
% all_actlon(jcount,1,:)=act_lon(:,:);
% all_actlat(jcount,1,:)=act_lat(:,:);

all_templon(jcount,1:totalicount,:)=St_templon(:,:);
all_templat(jcount,1:totalicount,:)=St_templat(:,:);
%all_temph(jcount,1:totalicount,:)=St_tempH(:,:);
all_tempt(jcount,1:totalicount,:)=St_tempT(:,:);
all_tempah(jcount,1:totalicount,:)=St_tempAH(:,:);
all_tempU(jcount,1:totalicount,:)=St_tempU(:,:);
all_tempV(jcount,1:totalicount,:)=St_tempV(:,:);

hatch_templon(jcount,1,:)=all_templon(jcount,totalicount,:);
hatch_templat(jcount,1,:)=all_templat(jcount,totalicount,:);
%hatch_temph(jcount,1,:)=all_temph(jcount,totalicount,:);
hatch_tempt(jcount,1,:)=all_tempt(jcount,totalicount,:);
hatch_tempah(jcount,1,:)=all_tempah(jcount,totalicount,:);
hatch_tempU(jcount,1,:)=all_tempU(jcount,totalicount,:);
hatch_tempV(jcount,1,:)=all_tempV(jcount,totalicount,:);
% all_dis(jcount)=dis_drift;
   


end        

save("2019_spawning_locations_3m_k_0.1_100_boundary2_10.mat",'hatch_tempU','hatch_tempV','hatch_tempah','hatch_templat','hatch_templon','hatch_tempt');
save("2019_spawning_locations_3m_k_0.1_100_trajectory_boundary2_10.mat",'all_tempU','all_tempV','all_tempah','all_templat','all_templon','all_tempt');
