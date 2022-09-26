clear all
close all
format long g 
fclose all

filenameAT='OSV4_AT.txt';
filenameET='OSV4_ET.txt';
filenameGW='OSV4_GW.txt';

rc=4*2.54/100/2; % Casing radius
Laq=9; % Aquifer thickness or screen length in case of partial penetretion well
gam=sqrt(0.5)/(Laq/rc); % in case of partial penetretion well

partial=0; % in case of partial penetretion well : partial =1;

nbit=100000 % iteration number of semi-random Monte Carlo exploration using Hvorsel model
minfit=1  % Maximum fit to plot the results (K, fit) in fig 1
plotfit=0.597; % Maximum fit to plot the impulses responses with good fit in fig 6001
imin=-7; %log10(K) range : minimum  10^(-5)
imax=-4; % %log10(K) range : max  10^(-3)

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%  open data  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirt=chdir();

fib=fopen(filenameAT,'r'); i=1;
if fib>0,     ligne = fgetl(fib);
    while length(ligne)>1                
          [nul] = sscanf(ligne,'%d/%d/%d %d:%d:%d,%f');
          B(i,1)=datenum([nul(3), nul(2),nul(1), nul(4), nul(5), nul(6)]);        
          B(i,2)=nul(7)/100;
          ligne = fgetl(fib);i=i+1;
    end, fclose(fib);    
end

fib=fopen(filenameET,'r'); i=1;
if fib>0,     ligne = fgetl(fib);
    while length(ligne)>1                
          [nul] = sscanf(ligne,'%d/%d/%d %d:%d:%d,%f');
          E(i,1)=datenum([nul(3), nul(2),nul(1), nul(4), nul(5), nul(6)]);        
          E(i,2)=nul(7);
          ligne = fgetl(fib);i=i+1;
    end, fclose(fib);    
end

fib=fopen(filenameGW,'r'); i=1;
if fib>0,    ligne = fgetl(fib);
    while length(ligne)>1                
          [nul] = sscanf(ligne,'%d/%d/%d %d:%d:%d,%f');
          D(i,1)=datenum([nul(3), nul(2),nul(1), nul(4), nul(5), nul(6)]);        
          D(i,2)=nul(7)/100;
          ligne = fgetl(fib);i=i+1;
    end,   fclose(fib);    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Filters  %%%%%%%%%%
dfilt1 = designfilt('highpassiir','StopbandFrequency',1, 'PassbandFrequency',1.5, 'StopbandAttenuation',55, 'PassbandRipple',4, 'DesignMethod','butter', 'MatchExactly','stopband','SampleRate',72)   ;            % Sample rate
%fvtool(dfilt1)
dfilt = designfilt('highpassiir','StopbandFrequency',0.2, 'PassbandFrequency',0.5, 'StopbandAttenuation',55, 'PassbandRipple',4, 'DesignMethod','butter', 'MatchExactly','stopband','SampleRate',72)  ;
%%%%%%%%%%%%%%%%%%%%%1 
 

yd1 = filtfilt(dfilt,D(:,2)); 
yb1 = filtfilt(dfilt,B(:,2)); 
ye1 = filtfilt(dfilt,E(:,2));  
 
 
F=[0.929536 0.997262 1 1.002738 1.895982 1.932274 2 2.005476 0.893234 ];
W=2*pi*F;
Wt=W;

nbc=length(W);           
F=F(1:nbc);                
W=W(1:nbc); 

Fs = (24*round(1/((D(2,1)-D(1,1))*24))); 



%%%%%%%%%%%%%%%%%%%%%%  HALS   %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%   ET  %%%%%%%%%%%%%%%%%%%%
tt=E(:,1); 
dd=ye1;

W=Wt([1 2 4:9]); nbc=length(W);

M=zeros(nbc*2,nbc*2);  Y=zeros(nbc*2,1);         
for i=1:2:nbc*2
Y(i)=sum(dd.*sin(W(floor((i+1)/2)).*tt));  
Y(i+1)=sum(dd.*cos(W(floor((i+1)/2)).*tt));
    for j=1:2:nbc*2
    M(i,j:j+1)=[sum(sin(W(floor((i+1)/2)).*tt).*sin(W(floor((j+1)/2)).*tt)),   sum(sin(W(floor((i+1)/2)).*tt).*cos(W(floor((j+1)/2)).*tt))]; 
    M(i+1,j:j+1)=[sum(cos(W(floor((i+1)/2)).*tt).*sin(W(floor((j+1)/2)).*tt)),   sum(cos(W(floor((i+1)/2)).*tt).*cos(W(floor((j+1)/2)).*tt))] ;  
    end
end

R=inv(M)*Y; ddd=zeros(length(dd),1);
for i=1:2:nbc*2,
    ddd=ddd+R(i).*sin(W(floor((i+1)/2)).*tt) + R(i+1).*cos(W(floor((i+1)/2)).*tt); %  +   R(3).*sin(Wm2.*tt) + R(4).*cos(Wm2.*tt)  +   R(5).*sin(Ws2.*tt) + R(6).*cos(Ws2.*tt)  +   R(7).*sin(Wo1.*tt) + R(8).*cos(Wo1.*tt)  +   R(9).*sin(Wn2.*tt) + R(10).*cos(Wn2.*tt)
end
disp('fit SourceET')
fit=sqrt(sum((dd-ddd).*(dd-ddd))/length(dd))/std(dd)

R1=R;
for i=1:2:length(R), 
    teta1t(floor(i/2)+1)=(atan(R(i+1)/R(i))); % teta2(floor(i/2)+1)=(atan(R(i)/R(1+i)))/W(floor(i/2)+1); 
    Ampl(floor(i/2)+1)=sqrt(R(i)*R(i)+R(i+1)*R(i+1));
    teta1s(floor(i/2)+1)=asin(R(i+1)/Ampl(floor(i/2)+1));
    teta1c(floor(i/2)+1)=acos(R(i)/Ampl(floor(i/2)+1)); 
    if R(i)<0 & R(i+1)<0, tetagood1(floor(i/2)+1)=2*pi-teta1c(floor(i/2)+1); end
    if R(i)<0 & R(i+1)>0, tetagood1(floor(i/2)+1)=teta1c(floor(i/2)+1); end
    if R(i)>0, tetagood1(floor(i/2)+1)=teta1s(floor(i/2)+1); end        
end

ddd1=zeros(length(dd),1);
for i=1:nbc,
    ddd1=ddd1+Ampl(i)*sin(W(i).*tt+ tetagood1(i)) ;
end

AmplET=[Ampl(1) Ampl(2) 0 Ampl(3) Ampl(4) Ampl(5) Ampl(6) Ampl(7:end)]; AmplET./max(AmplET)
tetaET=[tetagood1(1) tetagood1(2) -2*pi tetagood1(3) tetagood1(4) tetagood1(5) tetagood1(6) tetagood1(7:end)];
W=Wt; nbc=length(W);

%%%%%%%%%   BARO  %%%%%%%%%%%%%%%%%%%%
tt=B(:,1); 
dd=yb1;

W=Wt([ 2 3 4 5 6 7 8])% 9:length(Wt)]); 
nbc=length(W)

M=zeros(nbc*2,nbc*2);  Y=zeros(nbc*2,1);         
for i=1:2:nbc*2
Y(i)=sum(dd.*sin(W(floor((i+1)/2)).*tt));  
Y(i+1)=sum(dd.*cos(W(floor((i+1)/2)).*tt));
    for j=1:2:nbc*2
    M(i,j:j+1)=[sum(sin(W(floor((i+1)/2)).*tt).*sin(W(floor((j+1)/2)).*tt)),   sum(sin(W(floor((i+1)/2)).*tt).*cos(W(floor((j+1)/2)).*tt))]; 
    M(i+1,j:j+1)=[sum(cos(W(floor((i+1)/2)).*tt).*sin(W(floor((j+1)/2)).*tt)),   sum(cos(W(floor((i+1)/2)).*tt).*cos(W(floor((j+1)/2)).*tt))] ;  
    end
end

R=inv(M)*Y; ddd=zeros(length(dd),1);
for i=1:2:nbc*2,
    ddd=ddd+R(i).*sin(W(floor((i+1)/2)).*tt) + R(i+1).*cos(W(floor((i+1)/2)).*tt); %  +   R(3).*sin(Wm2.*tt) + R(4).*cos(Wm2.*tt)  +   R(5).*sin(Ws2.*tt) + R(6).*cos(Ws2.*tt)  +   R(7).*sin(Wo1.*tt) + R(8).*cos(Wo1.*tt)  +   R(9).*sin(Wn2.*tt) + R(10).*cos(Wn2.*tt)
end
disp('fit SourceAT')
fit=sqrt(sum((dd-ddd).*(dd-ddd))/length(dd))/std(dd)

R1=R; Ampl=[]; tetagood1=[];
for i=1:2:length(R), 
    teta1t(floor(i/2)+1)=(atan(R(i+1)/R(i))); % teta2(floor(i/2)+1)=(atan(R(i)/R(1+i)))/W(floor(i/2)+1); 
    Ampl(floor(i/2)+1)=sqrt(R(i)*R(i)+R(i+1)*R(i+1));
    teta1s(floor(i/2)+1)=asin(R(i+1)/Ampl(floor(i/2)+1));
    teta1c(floor(i/2)+1)=acos(R(i)/Ampl(floor(i/2)+1)); 
    if R(i)<0 & R(i+1)<0, tetagood1(floor(i/2)+1)=2*pi-teta1c(floor(i/2)+1); end
    if R(i)<0 & R(i+1)>0, tetagood1(floor(i/2)+1)=teta1c(floor(i/2)+1); end
    if R(i)>0, tetagood1(floor(i/2)+1)=teta1s(floor(i/2)+1); end        
end

ddd1=zeros(length(dd),1);
for i=1:nbc,
    ddd1=ddd1+Ampl(i)*sin(W(i).*tt+ tetagood1(i)) ;
end

%W=Wt([3 5 6 7 8 9:length(Wt)]);
AmplAT=[ 0 Ampl(1) Ampl(2) Ampl(3) Ampl(4) Ampl(5) Ampl(6:length(Ampl)) 0];  AmplAT./max(AmplAT)
tetaAT=[ 0 tetagood1(1) tetagood1(2) tetagood1(3) tetagood1(4) tetagood1(5)  tetagood1(6:length(Ampl)) 0];
W=Wt; nbc=length(W);


%%%%%%%%%   GW  %%%%%%%%%%%%%%%%%%%%
tt=D(:,1); 
dd=yd1;

W=Wt(1:9); nbc=length(W)

M=zeros(nbc*2,nbc*2);  Y=zeros(nbc*2,1);         
for i=1:2:nbc*2
Y(i)=sum(dd.*sin(W(floor((i+1)/2)).*tt));  
Y(i+1)=sum(dd.*cos(W(floor((i+1)/2)).*tt));
    for j=1:2:nbc*2
    M(i,j:j+1)=[sum(sin(W(floor((i+1)/2)).*tt).*sin(W(floor((j+1)/2)).*tt)),   sum(sin(W(floor((i+1)/2)).*tt).*cos(W(floor((j+1)/2)).*tt))]; 
    M(i+1,j:j+1)=[sum(cos(W(floor((i+1)/2)).*tt).*sin(W(floor((j+1)/2)).*tt)),   sum(cos(W(floor((i+1)/2)).*tt).*cos(W(floor((j+1)/2)).*tt))] ;  
    end
end

R=inv(M)*Y; ddd=zeros(length(dd),1);
for i=1:2:nbc*2,
    ddd=ddd+R(i).*sin(W(floor((i+1)/2)).*tt) + R(i+1).*cos(W(floor((i+1)/2)).*tt); %  +   R(3).*sin(Wm2.*tt) + R(4).*cos(Wm2.*tt)  +   R(5).*sin(Ws2.*tt) + R(6).*cos(Ws2.*tt)  +   R(7).*sin(Wo1.*tt) + R(8).*cos(Wo1.*tt)  +   R(9).*sin(Wn2.*tt) + R(10).*cos(Wn2.*tt)
end
disp('fit GW')
fit=sqrt(sum((dd-ddd).*(dd-ddd))/length(dd))/std(dd)

R1=R; Ampl=[]; tetagood1=[];
for i=1:2:length(R), 
    teta1t(floor(i/2)+1)=(atan(R(i+1)/R(i))); % teta2(floor(i/2)+1)=(atan(R(i)/R(1+i)))/W(floor(i/2)+1); 
    Ampl(floor(i/2)+1)=sqrt(R(i)*R(i)+R(i+1)*R(i+1));
    teta1s(floor(i/2)+1)=asin(R(i+1)/Ampl(floor(i/2)+1));
    teta1c(floor(i/2)+1)=acos(R(i)/Ampl(floor(i/2)+1)); 
    if R(i)<0 & R(i+1)<0, tetagood1(floor(i/2)+1)=2*pi-teta1c(floor(i/2)+1); end
    if R(i)<0 & R(i+1)>0, tetagood1(floor(i/2)+1)=teta1c(floor(i/2)+1); end
    if R(i)>0, tetagood1(floor(i/2)+1)=teta1s(floor(i/2)+1); end        
end

ddd1=zeros(length(dd),1);
for i=1:nbc,
    ddd1=ddd1+Ampl(i)*sin(W(i).*tt+ tetagood1(i)) ;
end

AmplGW=Ampl;
tetaGW=tetagood1;
%W=Wt; nbc=length(W);


%%%%%%%%%%%%%%%%%%%%%%%%  Disantgling %%%%%%%%%%%%
tic()
at=AmplGW(3)*sin(W(3)*tt+tetaGW(3))   +   AmplGW(7)*sin(W(7)*tt+tetaGW(7))  -   AmplET(7)/AmplET(6)*AmplGW(6)*sin(W(7)*tt+tetaGW(6)-tetaET(6)+tetaET(7));   +   AmplGW(5)*sin(W(5)*tt+tetaGW(5))  -   AmplET(5)/AmplET(6)*AmplGW(6)*sin(W(5)*tt+tetaGW(6)-tetaET(6)+tetaET(5));
%at=AmplGW(3)*sin(W(3)*tt+tetaGW(3))-AmplET(3)/AmplET(2)*AmplGW(2)*sin(W(3)*tt+tetaGW(2)-tetaET(2)+tetaET(3))   +   AmplGW(1)*sin(W(1)*tt+tetaGW(1))-AmplET(1)/AmplET(4)*AmplGW(4)*sin(W(1)*tt+tetaGW(4)-tetaET(4)+tetaET(1));

dd=at;

M=zeros(nbc*2,nbc*2);
Y=zeros(nbc*2,1);

for i=1:2:nbc*2
Y(i)=sum(dd.*sin(W(floor((i+1)/2)).*tt));  
Y(i+1)=sum(dd.*cos(W(floor((i+1)/2)).*tt));
    for j=1:2:nbc*2
    M(i,j:j+1)=[sum(sin(W(floor((i+1)/2)).*tt).*sin(W(floor((j+1)/2)).*tt)),   sum(sin(W(floor((i+1)/2)).*tt).*cos(W(floor((j+1)/2)).*tt))]; 
    M(i+1,j:j+1)=[sum(cos(W(floor((i+1)/2)).*tt).*sin(W(floor((j+1)/2)).*tt)),   sum(cos(W(floor((i+1)/2)).*tt).*cos(W(floor((j+1)/2)).*tt))] ;  
    end
end

R=inv(M)*Y; ddd=zeros(length(dd),1);
for i=1:2:nbc*2,
    ddd=ddd+R(i).*sin(W(floor((i+1)/2)).*tt) + R(i+1).*cos(W(floor((i+1)/2)).*tt); %  +   R(3).*sin(Wm2.*tt) + R(4).*cos(Wm2.*tt)  +   R(5).*sin(Ws2.*tt) + R(6).*cos(Ws2.*tt)  +   R(7).*sin(Wo1.*tt) + R(8).*cos(Wo1.*tt)  +   R(9).*sin(Wn2.*tt) + R(10).*cos(Wn2.*tt)
end
%disp('fit GWL S2 (AT)')
fit=sqrt(sum((dd-ddd).*(dd-ddd))/length(dd));

R2=R;
for i=1:2:length(R), 
    teta2t(floor(i/2)+1)=(atan(R(i+1)/R(i))); % teta2(floor(i/2)+1)=(atan(R(i)/R(1+i)))/W(floor(i/2)+1); 
    Ampl3(floor(i/2)+1)=sqrt(R(i)*R(i)+R(i+1)*R(i+1));
    teta2s(floor(i/2)+1)=asin(R(i+1)/Ampl3(floor(i/2)+1));
    teta2c(floor(i/2)+1)=acos(R(i)/Ampl3(floor(i/2)+1)); 
    if R(i)<0 & R(i+1)<0, tetagood3(floor(i/2)+1)=2*pi-teta2c(floor(i/2)+1); end
    if R(i)<0 & R(i+1)>0, tetagood3(floor(i/2)+1)=teta2c(floor(i/2)+1); end
    if R(i)>0, tetagood3(floor(i/2)+1)=teta2s(floor(i/2)+1); end      
end

ddd1=zeros(length(dd),1); ddd2=zeros(length(dd),1);
for i=1:nbc,
    ddd1=ddd1+Ampl3(i)*sin(W(i).*tt+ tetagood3(i)) ;
end

disp('dephasage GW-ET (M2) °')
(tetaGW(6)-tetaET(6))*180/pi

disp('dephasage GWat-Baro (S1) °')
(tetaGW(3)-tetaAT(3))*180/pi

disp('BARO S1')
AmplGW(3)/AmplAT(3)'

disp('dephasage GWat-Baro (S2) °')
dephaGWATs2=(tetagood3(7)-tetaAT(7))*180/pi

disp('BARRROOOOO formula RAU')
rau=Ampl3'./AmplAT'; 
Ampl3
if max(rau([1,2,3,4,5,6,8]))>0.01*rau(7), rau, disp('dephasage GWat-Baro (S1) °'), (tetagood3(1)-tetaAT(1))*180/pi

else, rau(7),
end
toc()
%%%%%%  Disantangling general en assumant que les semis-diurnaux répondent
%%%%%%  pareil en phase & ampl, mais bien sur different pour ET et AT
yd1 = filtfilt(dfilt1,D(:,2)); 

fun = @(x)  sqrt(sum( (-yd1 + x(1)*(AmplAT(5)*sin(W(5).*tt+tetaAT(5)+x(2)) + AmplAT(6)*sin(W(6).*tt+tetaAT(6)+x(2)) + AmplAT(7)*sin(W(7).*tt+tetaAT(7)+x(2)) + AmplAT(8)*sin(W(8).*tt+tetaAT(8)+x(2)))   +   x(3)*(AmplET(5)*sin(W(5).*tt+tetaET(5)+x(4)) + AmplET(6)*sin(W(6).*tt+tetaET(6)+x(4)) + AmplET(7)*sin(W(7).*tt+tetaET(7)+x(4)) + AmplET(8)*sin(W(8).*tt+tetaET(8)+x(4))) ).^2)/length(tt))
x0 = [0.58 200*pi/180 0.00028 10*pi/180];
%[x, fval] = fminsearch(fun,x0)
[x, fval] = fminsearch(fun,[0.8 0 0.0001 0 ])
if x(1)<0,  [x, fval] = fminsearch(fun,[-x(1) x(2)+pi x(3) x(4) ]), end

disp('Rau vs solver method')
disp('BE'), [rau(7) x(1)], be=x(1);
disp('BE should be betweent [0 1] ')
disp('Depha AT'), [ dephaGWATs2 x(2)*180/pi-360], dat=x(2)*180/pi-360;
disp('WS'), [AmplGW(6)/AmplET(6) x(3)], ws=x(3);
disp('Depha ET'), [ (tetaGW(6)-tetaET(6))*180/pi   x(4)*180/pi], det=x(4)*180/pi;
x2=x; x00=x;

stepd=round((D(2)-D(1))*24*60);  
disp('data Time step (min)      AT S2 depha (min)    ET M2 depha (min) ')
[stepd, (dat+180)/(360)*12*60 det/(360)*24/F(6)*60]
disp('AT S2 depha (min) should be cromprised between [-1 1]*data Time step (min) ')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INVERSION %

disp('Time step data (min)')
stepd=round((D(2)-D(1))*24*60)

yd1 = filtfilt(dfilt1,D(:,2)); 
yb1 = filtfilt(dfilt1,B(:,2)); 
ye1 = filtfilt(dfilt1,E(:,2));

   
    stepdd=stepd; % in minutes
    
    dtmp=D(1,1):stepdd/(24*60):D(end,1); tt=dtmp'; 
    tmp = interp1(D(:,1),yd1,dtmp,'pchip'); yd1c=tmp';
    tmp = interp1(D(:,1),yb1,dtmp,'pchip'); yb1c=tmp';
    tmp = interp1(D(:,1),ye1,dtmp,'pchip'); ye1c=tmp';
    
    %%%%% Defining response sizes
    dH=yd1c(1:1:end);
       
    iH=[]; iHi=0;
    while stepdd<6*60
        iH=[iH, stepdd]; iHi=[iHi, stepdd/stepd];
        stepdd=2*stepdd;
    end
      
%     tmp=((tetagood3(7)-tetaAT(7))*180/pi+180); if tmp<0, tmp=-tmp, elseif tmp>360, tmp=tmp-360; else tmp=tmp; end, tmp
%     round(tmp/(360*F(7)*stepdd/(24*60)))
%     n=max([round(round(tmp/(360*F(7)*stepdd/(24*60)))*2), 6])

    n=length(iHi);
    
    n2=max([(round(det/(360*F(7)*stepdd/(24*60)))+1)*2 3]);  % if n2>n,n=n2;end
        
    %%%% La matrice %%%

    %dH=dH(n+1:1:end-n2); tt=tt(n+1:1:end-n2);
    dH=dH(iHi(end)+1:1:end-n2); tt=tt(iHi(end)+1:1:end-n2);
    
    
    M=[];
    for i=1:length(dH)
        for j=1:n
            M(i,0*n+j)=yb1c(iHi(n)+i-iHi(j));   
        end
        for j=n+1:n+n2*2+1
            M(i,j)=ye1c(iHi(n)+i +n2-(j-(n+1)));       
        end
    end
    
    
 
    %%%%%%% MC-Hvorslev
    fitkk=[]; kk=0; fitmin=100; mu=[]; mu1=zeros(1,n);
              
        for i=1:nbit
           mu1=[]; nu1=[];

        %%%% model Horslev direct
         m= 1/(-rc^2*log(200)/(2*Laq)/(stepd*60)) * 10^(imax+(imin-imax)*rand());%  *Ki(ik)* (1+8*rand())^(sign(tmp1));       %-10^(-3+5*rand());
         bb=4*(rand()-0.5); %   if GW response start just before or just after the AT slug (due to low sampling for example)
         bb=0;
         
         if bb<=0, mu1(1)=0; for im=2:length(iHi), mu1(im)=(exp(m*(iHi(im)+bb))-exp(m*(iHi(im-1)+bb)))*be; end, end  % bb<0 demarre + tard
         if bb>0, mu1(1)=-(1-exp(m*bb))*be; for im=2:length(iHi), mu1(im)=(exp(m*(iHi(im)+bb))-exp(m*(iHi(im-1)+bb)))*be; end, end
             
%          if bb>=0,mu1(1:floor(bb)+1)=0;   mu1(floor(bb)+2:n)=(exp(m*(iHi([2:n-floor(bb)])-(bb)))-exp(m*(iHi([1:n-1-floor(bb)])-(bb))))*rau(7);        % demarre + tard
%          elseif bb<0,  mu1(1:n)=(exp(m*(iHi([0:n-1])-(bb)))-exp(m*(iHi([-1:n-2])-(bb))))*rau(7);     % demarre + tot
%          end
         
        %%%% model ET
           %npic=floor(rand()*(2*n2))+1;
           npic=n2+1-round(det/(360*F(7)*stepdd/(24*60)));  % Localized Apex according to ET-GW phase shift
           hpic=rand()*ws*1+0.3*ws;
           maxtmp=hpic; nu1(npic)=hpic;
           if npic>1,for j=npic-1:-1:1
               nu1(j)=rand()*maxtmp; maxtmp=nu1(j);
           end,end
           maxtmp=hpic;
           for j=npic+1:1:2*n2+1
               nu1(j)=rand()*maxtmp; maxtmp=nu1(j);
           end    

        mu=[mu1 nu1]'; mu2=nu1;  
        fit=sqrt(sum((dH-M*mu).*(dH-M*mu))/length(dH))/std(dH);
        fitkkk(i)=fit;

        mu1=mu(1:n); mu2=mu(n+1:end); 
       
  
           if fit<fitmin, fitmin=fit; fitminK=fit;  mug3=mu; mg3=m; bbg=bb; ydcg3=M*mug3;  end 
           
           if fit<minfit, 
               kk=kk+1; fitkk(kk)=fit;  Kfit(kk)= -rc^2*log(200)/(2*Laq)*m/(stepd*60); 
                Kfitparcial(kk)=-rc^2*log(1/(2*gam)+sqrt(1+1/(2*gam)^2))/(2*Laq)*1/(stepd*60)*m;
           end
                      
           if fit<plotfit, figure(6000), subplot(211), hold on, plot([0 iH],mu1),    subplot(212), hold on, plot([-n2:n2]/(60/stepd),mu2), end 
                                   
        end


fitminK
disp('HVORSLEV')
disp(['Kfit   ',num2str(-rc^2*log(200)/(2*Laq)*mg3/(stepd*60))])

mu=mug3;
dHat=M(:,1:n)*mu(1:n); disp('%AT in GW'),  mean(dHat./dH)
dHet=M(:,n+1:end)*mu(n+1:end); disp('%ET in GW'), mean(dHet./dH)

mu1=mu(1:n); mu2=mu(n+1:end); 
Smu1=zeros(n,1); for i=1:n, Smu1(i)=sum(mu1(1:i));end,   Smu2=zeros(2*n2+1,1); for i=1:2*n2+1, Smu2(i)=sum(mu2(1:i));end
figure(7001), subplot(2,1,1), hold on, plot([0 iH],Smu1,'b','DisplayName',['AT-S2-BE-Hvorslev K= ',num2str(-rc^2*log(200)/(2*Laq)*mg3/(stepd*60)),' m/s']),
xlabel(['Time']),ylabel('m/m'),title('semi-diurnal AT cumulative response'), legend
subplot(2,1,2), hold on, plot([-n2:n2]/(60/stepd),Smu2,'b','DisplayName','Inv MC-Hvorslev'),title('semi-diurnal ET cumulative response'), xlabel(['Time']),ylabel('m/nstrain'),

figure(6000), subplot(2,1,1), hold on, plot([0 iH],mu1,'b','DisplayName',['MC Hvorslev K= ',num2str(-rc^2*log(200)/(2*Laq)*mg3/(stepd*60)),' m/s']),
xlabel(['Time']),ylabel('m/m'),title('semi-diurnal AT impulse response'), , legend
subplot(2,1,2), hold on, plot([-n2:n2]/(60/stepd),mu2,'b','DisplayName','MC-Hvorslev'),title('semi-diurnal ET impulse response'), xlabel(['Time']),ylabel('m/nstrain'),
axis([-n2/(60/stepd) n2/(60/stepd) -0.1*ws 1.3*ws ])

%%% Fit of WS,BE and phases shifts
ydc=ydcg3; %tt=D(1,1):stepdd/(24*60):D(end,1);  tt=tt(n+1:1:end-n2); tt=tt';
fun = @(x)  sqrt(sum( (-ydc + x(1)*(AmplAT(5)*sin(W(5).*tt+tetaAT(5)+x(2)) + AmplAT(6)*sin(W(6).*tt+tetaAT(6)+x(2)) + AmplAT(7)*sin(W(7).*tt+tetaAT(7)+x(2)) + AmplAT(8)*sin(W(8).*tt+tetaAT(8)+x(2)))   +   x(3)*(AmplET(5)*sin(W(5).*tt+tetaET(5)+x(4)) + AmplET(6)*sin(W(6).*tt+tetaET(6)+x(4)) + AmplET(7)*sin(W(7).*tt+tetaET(7)+x(4)) + AmplET(8)*sin(W(8).*tt+tetaET(8)+x(4))) ).^2)/length(tt))
[x, fval] = fminsearch(fun,[x00(1) x00(2) x00(3) x00(4) ]);
if x(1)<0,  [x, fval] = fminsearch(fun,[-x(1) x(2)+pi x(3) x(4) ]); end
fitf2=sqrt(((x(1)-x00(1))/x00(1))^2+((x(3)-x00(3))/x00(3))^2+((x(2)-x00(2))/x00(2))^2+((x(4)-x00(4))/x00(4))^2/4);
disp(['fit MC Hvorslev(temp, freq):   ',num2str(fit),'   ',num2str(fitf2)])
fitminf2=fitf2; x1=x;
'raw,  fitH,  fitf2'
[x00; x1;x ]'

if partial==0
figure(1)
semilogx(Kfit,fitkk,'.k','DisplayName','Inv. MC-Hvorslev')
title('Results of BE-Hvorslev model for S2-AT slug test'), xlabel('K (m/s)'), ylabel('GWL fit RMSE (m)')
else
figure(2)
semilogx(Kfitparcial,fitkk,'.k','DisplayName','Inv. MC-Hvorslevpartialpenetration')
title('Results of BE-Hvorslev model for S2-AT slug test (partial penetrating well)'), xlabel('K (m/s)'), ylabel('GWL fit RMSE (m)')
end

