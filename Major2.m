%% Loading ecg signal

ecg_signal=load('C:\Users\priya\Downloads\100m.mat', 'val');
new_signal=ecg_signal.val(1,:);

t=0:1:3599;

figure(1)
plot(t,new_signal,'k-')
xlabel ('Time')
ylabel ('Amplitude')
legend ('ECG Signal')

signal=zeros(1,100);

for i=1:100
    signal(1,i)=new_signal(1,i)./200;
end
figure(2)
plot(signal,'k-')
xlabel ('Time')
ylabel ('Amplitude')
legend ('ECG Signal(1000 samples)')
axis([0 100 -0.6 1])


%% Adding noise to the signal

sampling_frequency = 1000;
mains_coeff = 0.3;   
time_step = 1/sampling_frequency;
max_time = 0.100;   
t = time_step:time_step:max_time; 

mains_signal = cos(2.*pi.*60.*t);

figure(3)
plot(mains_signal,'k-')
xlabel ('Time')
ylabel ('Amplitude')
legend ('Noise Signal')

noise_signal=mains_coeff*mains_signal;
dirty_signal=signal+noise_signal;

figure(4)
plot(dirty_signal,'k-')
xlabel ('Time')
ylabel ('Amplitude')
legend ('Noise Added Signal')


%% Bandpass filter

d = designfilt('bandpassfir', 'FilterOrder', 3, 'CutoffFrequency1', 5, 'CutoffFrequency2', 22, 'SampleRate', 360);
filtered_signal=filter(d,dirty_signal);
t=1:100;
figure(5)
plot(t,filtered_signal,'k-')
xlabel ('Time')
ylabel ('Amplitude')
legend ('Signal After it is passed through BPF')


%% Moth Flame Algorithm

T=500;
lb=-1;
ub=1;
n=30;
d=3;
M=unifrnd(lb,ub,n,d);
iteration=1;
while iteration<T+1
    %flameno=round(n-iteration*((n-1)/T));
    for i=1:n
        forub=M(i,:)>ub;
        forlb=M(i,:)<lb;
        M(i,:)=(M(i,:).*(~(forub+forlb)))+ub.*forub+lb.*forlb;
         OM(1,i)=DFOD(M(i,:));
    end
     if iteration==1
        [OMS I]=sort(OM);
        MS=M(I,:);
        BFS=MS;
        OBFS=OMS;
     else
         
        M2=[lastM;BFS];
        OM2=[lastOM OBFS];
        
        [OM2S I]=sort(OM2);
        M2S=M2(I,:);
        
        OMS=OM2S(1:n);
        MS=M2S(1:n,:);
        BFS=MS;
        OBFS=OMS;
    end
 
    OBF=OMS(1);
    BFP=MS(1,:);
      
    lastM=M;
    lastOM=OM;
    
%      a=-1+iteration*((-1)/T);
     
     for i=1:size(M,1)
        
        for z=1:size(M,2)
                M(randi(size(M,1)*size(M,2),1))=M(randi(size(M,1)*size(M,2),1))+randn(1)./4;
%             if i<=flameno
%                 Dis=abs(MS(i,z)-M(i,z));
%                 b=1;
%                 t=(a-1)*rand+1;
%                 M(i,z)=Dis*exp(b.*t).*cos(t.*2*pi)+MS(i,z);
%             end
%             if i>flameno
%                 Dis=abs(MS(i,z)-M(i,z));
%                 b=1;
%                 t=(a-1)*rand+1;
%                 M(i,z)=Dis*exp(b.*t).*cos(t.*2*pi)+MS(flameno,z);
%             end
        end
     end
     iteration=iteration+1;
end


%% Filtering using DFOD

arr=M(1,:)
arr1=arr(1,1:2);
arr2=arr(1,3);
D=filter(arr1,arr2,filtered_signal);

figure(6)
t=1:100;
plot(t,D,'k-')
xlabel ('Time')
ylabel ('Amplitude')
legend ('Output of DFOD')


%% Amplitude Normaliser

a=zeros(1,100);
for i=1:100
    
    a(i)=D(i)./max(abs(D));
    
end

figure(7)
t=1:100;
plot(t,a,'k-')
xlabel ('Time')
ylabel ('Amplitude')
legend ('Normalised Signal')


%% Smooth Waveform Generator

S=zeros(1,100);
for i=1:100
    
    S(i)=-(a(i).^2.*log((a(i)).^2));
    
end

% figure(8)
t=1:100;
% plot(t,S,'k-')
% xlabel ('Time')
% ylabel ('Amplitude')
% legend ('Smooth Waveform')

h=ones(1,31)/31;
delay=15;
signal1=conv(S,h);
signal1=signal1(15+[1:100]);
signal1=signal1/max(abs(signal1));
% figure(9)
% plot(t,signal1,'k-')
% xlabel ('Time')
% ylabel ('Amplitude')
% legend ('Moving Window Integration')



%% R-peak finding logic using thresholding

vth=.85*max(D);
R_peak=0;

for j=1:1:100
    if D(1,j)>vth
        
        if R_peak>=D(1,j)
            R_peak=R_peak;
            
        else
            R_peak=D(1,j);
        end
         
    end
end

Identified_R_peak=R_peak


%% RLE Coding:

x=find(D==R_peak)
Binary_Matrix=zeros(1,100);
Binary_Matrix(x)=1;

figure(10)
plot(Binary_Matrix,'k-')
xlabel ('Time')
ylabel ('Amplitude')
legend ('Binary Matrix For RLE')

a=Binary_Matrix;
g=size(a);
h=g(2);
k=0;
j=zeros(1,h);
for i=1:h
    
    if(a(i)==1 & i==1)
         if(a(i+1)==0)
             k=k+1;
                j(i)=k.*10+1;
                   
               
        else
            k=k+1;
               
         end
         
    elseif(a(i)==1 & i==g(2))
        if(a(i-1)==1)
            k=k+1;
            j(i)=k.*10+1;
        else
            k=0;
            k=k+1;
            j(i)=k.*10+1;
        end
            
    
        elseif(a(i)==1 & i~=g(2))
             if(a(i-1)==1)
                if(a(i+1)==0)
                    k=k+1;
                        j(i)=k.*10+1;
                        
                else
                    k=k+1;
                end
                
            else
                k=0;
                
                    if(a(i+1)==0)
                    k=k+1;
                        j(i)=k.*10+1;
                        
                    else
                        k=k+1;
                        
                    end
    
                end
                
    else
        if(a(i)==0 & i==1)
            if (a(i+1)==1)
                k=k+1;
                    j(i)=k.*10;
                    
            else
                k=k+1;
                    
            end
            
       elseif(a(i)==0 & i==g(2))
           if(a(i-1)==0)
              k=k+1;
              j(i)=k.*10;
          else
            k=0;
            k=k+1;
            j(i)=k.*10;
        end
            
        elseif(a(i)==0 & i~=g(2))
             if(a(i-1)==0)
                if(a(i+1)==1)
                    k=k+1;
                        j(i)=10.*k;
                        
                else
                    k=k+1;
                end
            else
                k=0;
                    if(a(i+1)==1)
                    k=k+1;
                        j(i)=k.*10;
                        
                    else
                        k=k+1;
                    end
            end
        end
    end
    
    
end

RLE=j;
z=1;
Result=zeros(1,2);
for i=1:100
    
    if RLE(i)~=0
        Result(z)=RLE(i);
        z=z+1;
    end
    
end

RLE_coded_data=Result
