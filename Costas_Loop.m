%===================================================%
%              Simulation Homework Wk.10            %
%                  -信卓1601班—卫天峻-               %
%                   -ID:U201613509-                 %
%                -----Costas Loop-----              %
%===================================================%
function Costas_Loop
%% Parameters
% LO = 10.0005MHz
% RF = 10.0000MHz
% K_VCO = 4000
% Sample Rate = 1.00GHz
% Sample Length = 10000
Sp_Rate = 1e9;                                  %Sample Rate
Sp_Len = 1e4;                                   %Sample Length
Time = 0:1/Sp_Rate:(Sp_Len/Sp_Rate)-1/Sp_Rate;  %Time Scale
LO = 10.0005e6;                                 %Local Oscillate Frequency
RF = 10e6;                                      %Radio Frequency
K_VCO = 4e3;                                    %Parameter for voltage-controlled oscillator 
Theta_Error = 0;                                %Voltage-controlled oscillator Theta Error
%% Input Signal
Message = 1.5+cos(2*pi*3e5.*Time) +...          %Message Signal
            0.5*cos(2*pi*2e5.*Time);
Carrier = cos(2*pi*RF.*Time);                   %Carrier Signal
DSB_SC = awgn(Message.*Carrier,35,'measured');  %35dB Gaussian White Noise
%% Filter
lowpass_Filt = designfilt('lowpassfir',...      %Lowpass Filter
                'PassbandFrequency',1e5, 'StopbandFrequency',5e6,...
                'PassbandRipple',0.5,'StopbandAttenuation',45,...
                'DesignMethod','kaiserwin','SampleRate',1e9);
loop_Filt = [1 (1/200)*ones(1,200)];            %Loop Filter
%% Initialize
VLO_I = zeros(1,length(Time));                  %VLO(Initial)   
VLO_Q = zeros(1,length(Time));                  %VLO(Quadrature)
VI = zeros(1,length(Time));                     %Vout(Initial)
VQ = zeros(1,length(Time));                     %Vout(Quadrature)
V_PD = zeros(1,length(Time));                   %V Baseband Detector
V_OP = zeros(1,length(Time));                   %V Operation
ThetaE_Rec = zeros(1,length(Time));          %Theta Error Record
%% Loop
disp('Costas Loop Processing...');
for Cnt = 1:length(Time)
    VLO_I(Cnt) = 2*cos(2*pi*(LO+Theta_Error)*Time(Cnt));
    VLO_Q(Cnt) = 2*sin(2*pi*(LO+Theta_Error)*Time(Cnt));
    VI(Cnt) = filter(lowpass_Filt,VLO_I(1:Cnt).*DSB_SC(1:Cnt))*[zeros(Cnt-1,1);1];
    VQ(Cnt) = filter(lowpass_Filt,VLO_Q(1:Cnt).*DSB_SC(1:Cnt))*[zeros(Cnt-1,1);1];
    V_PD(Cnt) = VI(Cnt)*VQ(Cnt);
    V_OP(Cnt) = filter(loop_Filt(2),loop_Filt(1),V_PD(1:Cnt))*[zeros(Cnt-1,1);1];
    Theta_Error = Theta_Error - K_VCO*V_OP(Cnt);
    ThetaE_Rec(Cnt) = Theta_Error;
end
%% Graph
Freq = (0:Sp_Len-1) .* (Sp_Rate/Sp_Len);    %Frequency Domain
figure('Color','white','Name','Costas Loop - InputSignal','NumberTitle','off');
plot(subplot(3,3,[1 2]),Message);   Graph_Set(gca,[0 Sp_Len],[0 4],'Time/ns','Volt/V','Message Signal');
plot(subplot(3,3,[4 5]),DSB_SC);    Graph_Set(gca,[0 Sp_Len],[-6 6],'Time/ns','Volt/V','DSB-SC Signal');
plot(subplot(3,3,[7 8]),Carrier);   Graph_Set(gca,[0 Sp_Len],[-2 2],'Time/ns','Volt/V','Carrier Signal'); 
stem(subplot(3,3,3),Freq,abs(fft(Message)),'Marker','.');    Graph_Set(gca,[0 10e5],[0 2e4],'Freq/Hz','Magnitude','Message Signal');
stem(subplot(3,3,6),Freq,abs(fft(DSB_SC)),'Marker','.');     Graph_Set(gca,[9e6 11e6],[0 2e4],'Freq/Hz','Magnitude','DSB-SC Signal');
stem(subplot(3,3,9),Freq,abs(fft(Carrier)),'Marker','.');    Graph_Set(gca,[9e6 11e6],[0 2e4],'Freq/Hz','Magnitude','Carrier Signal');
figure('Color','white','Name','Costas Loop - Process','NumberTitle','off');
plot(subplot(2,2,1),VLO_I);   Graph_Set(gca,[0 Sp_Len],[-3 3],'Time/ns','Volt/V','Local Oscillate Initial');
plot(subplot(2,2,3),VLO_Q);   Graph_Set(gca,[0 Sp_Len],[-3 3],'Time/ns','Volt/V','Local Oscillate Quadrature');
plot(subplot(2,2,2),1:Sp_Len,VI,1:Sp_Len,Message);   Graph_Set(gca,[0 Sp_Len],[0 4],'Time/ns','Volt/V','Output Initial Compared with Message');
plot(subplot(2,2,4),VQ);   Graph_Set(gca,[0 Sp_Len],[-0.15 0.15],'Time/ns','Volt/V','Output Quadrature');
figure('Color','white','Name','Costas Loop - Process & Result','NumberTitle','off');
plot(subplot(2,6,[1 4]),V_PD);      Graph_Set(gca,[0 Sp_Len],[-3e-2 3e-2],'Time/ns','Volt/V','Baseand Detector');
stem(subplot(2,6,[5 6]),Freq,abs(fft(V_OP)),'Marker','.');    Graph_Set(gca,[0 3e7],[0 0.3],'Freq/Hz','Magnitude','Baseand Detector');
plot(subplot(2,6,[7 10]),V_OP);     Graph_Set(gca,[0 Sp_Len],[-2e-4 2e-4],'Time/ns','Volt/V','Operation Volt');
plot(subplot(2,6,[11 12]),1:Sp_Len,ThetaE_Rec,1:Sp_Len,ones(1,Sp_Len)*(RF-LO));      Graph_Set(gca,[0 Sp_Len],[-1e3 0],'Time/ns','Volt/V','Theta Error');
disp('Final Theta Error:    '+string(ThetaE_Rec(end)));

function Graph_Set(ax,XLim,YLim,Xlabel,Ylabel,Title)
%% Graph Setting Helper 
ax.FontSize = 11;
ax.Title.String = Title;
ax.XLabel.String = Xlabel;
ax.YLabel.String = Ylabel;
ax.XLim = XLim;
ax.YLim = YLim;
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
