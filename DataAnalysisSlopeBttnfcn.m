close all;clear;
Num_bead = 0;
%open motor file
[FileName_motor,PathName_motor] = uigetfile('D:\Program files\Matlab\myCode\.txt','Select a txt file');
[fid,message] = fopen(FileName_motor,'rt');
if fid == -1
    disp (message);
else
    disp (fid);
end
DELIMITER = '\t';
HEADERLINE = 0;
Motor = importdata(FileName_motor,DELIMITER,HEADERLINE);
fclose(fid);
MagPos = [Motor(:,1) Motor(:,3)];
%open Zdet file
[FileName_XYZ,PathName_bead] = uigetfile('D:\Program files\Matlab\myCode\.txt','Select a txt file');
[fid,message] = fopen(FileName_XYZ,'rt');
if fid == -1
    disp (message);
else
    disp (fid);
end
XYZ = textscan(fid,'','delimiter','\t','HeaderLines',1,'TreatAsEmpty',{'1.#QNAN00'});
fclose(fid);
for i = 1:length(XYZ)
    if isnan(XYZ{i}(end))
        XYZ{i} = zeros(length(XYZ{i}),1);
    end
end
Smooth = 100;
Num_totalBead = fix((length(XYZ) - 5)/3);
Num_reference = Num_totalBead - 1;% 0 represents the last bead
framerate = 500;
% for i = 5:3:5 + 3*Num_totalBead
%     %XYZsub{i} = - XYZ{i} + sgolayfilt(XYZ{5+ 3*Num_reference},4,Smooth - 1);
%     XYZsub{i} = - XYZ{i};
% end
XYZsub = -XYZ{5 + Num_bead*3};
Frame2 = 32;
for i = 1:round(length(XYZsub)/Frame2) - 1
    XYZsubFrame2(i,1) = i;
    XYZsubFrame2(i,2) = mean(XYZsub(i*Frame2 - Frame2 + 1:i*Frame2));
end
XYZSmooth = sgolayfilt(XYZsub,4,Smooth - 1);
XYZsubDiff = diff(XYZsub);
for NumBreakStart = 1:length(XYZsubDiff)
    if abs(XYZsubDiff(NumBreakStart)) >= 1
        NumBreakStart = NumBreakStart - 1;
        break;
    end
end
for NumBreakEnd = length(XYZsubDiff):-1:1
    if abs(XYZsubDiff(NumBreakEnd)) >= 1
        NumBreakEnd = NumBreakEnd + 1;
        break;
    end
end
fig = figure('units','normalized','outerposition',[0 0.2 1 0.8]);
p = plot(XYZ{1}/framerate,XYZsub,'s','MarkerSize',3,'Color',[0 0.4470 0.7410],'hittest', 'off');hold on;
p = plot(XYZsubFrame2(:,1)*Frame2/framerate,XYZsubFrame2(:,2),'o-','MarkerSize',3,'Linewidth',1,'Color',[0.8500 0.3250 0.0980],'hittest', 'off');hold on;
%p = plot(XYZ{1}(1:end - 1)/framerate,XYZsubDiff,'Linestyle','-','Linewidth',1,'Color',[0.8500 0.3250 0.0980],'hittest', 'off');hold on;
%p = plot(XYZ{1}/framerate,XYZSmooth,'Linestyle','-','Linewidth',1,'Color',[0.8500 0.3250 0.0980],'hittest', 'off');hold on;
set(gca,'position',[0.05 0.11 0.9 0.815],'color',[1 1 1]);
grid on;
grid minor;
ax = gca;
ax.GridAlpha = 1;
ax.MinorGridAlpha = 1;
ax.GridColor = [0.7 0.7 0.7];
ax.MinorGridColor = [0.7 0.7 0.7];
title(strcat('DNA extension vs time-',FileName_XYZ(1:end - 4),'-bead-',num2str(Num_bead)));
xlabel('Time (s)');
ylabel('DNA extension (um)');
%ylim([min(min(XYZsub(1:NumBreakStart),min(XYZsub(NumBreakEnd:end)))) - 0.2,max(max(XYZsub(1:NumBreakStart)),max(XYZsub(NumBreakEnd:end))) + 0.2]);
ylim([mean(XYZsub) - 4, mean(XYZsub) + 3]);
xlim([min(XYZ{1}/framerate) - 50, max(XYZ{1}/framerate) + 50]);
if exist('Points.txt', 'file')
    delete Points.txt;
end
set(gca,'buttondownfcn',@mybttnfcn_);

%open the Points file, calculate the duration and amplitude of each segment
%and save them.
[FileName_Save,PathName_Save] = uigetfile('*.txt','Select a txt file');
[fid,message] = fopen(FileName_Save,'rt');
if fid == -1
    disp (message);
else
    disp (fid);
end
Points = textscan(fid,'','delimiter','\t','HeaderLines',0,'TreatAsEmpty',{'1.#QNAN00'});
fclose(fid);
k = 1;
for i = 1:length(Points{3}) - 1
    if Points{3}(i + 1) - Points{3}(i) > 0
        Save(k,1) = round(Points{1}(i)*framerate)/framerate;% x value in seconds of the first point
        Save(k,2) = XYZSmooth(round(Points{1}(i)*framerate));% y value in nanometer of the first point
        Save(k,3) = round(Points{1}(i + 1)*framerate)/framerate;% x value in seconds of the second point
        Save(k,4) = XYZSmooth(round(Points{1}(i + 1)*framerate));% y value in nanometer of the second point
        Save(k,5) = Save(k,3) - Save(k,1);% delta x in seconds
        Save(k,6) = Save(k,4) - Save(k,2);% delta y in nanometers
        Save(k,7) = Save(k,6)/Save(k,5);% slope (delta y over delta x)
        Save(k,8) = Points{3}(i);% segment sequence
        Save(k,9) = MagPos(round(framerate*Points{1}(i)),2);% magnets position
        fitDuration = (round(Points{1}(i)*framerate):round(Points{1}(i + 1)*framerate))';        
        Linearfit = polyfit(fitDuration/framerate,XYZsub(fitDuration),1);
        Save(k,10) = Linearfit(1);% slope 
        Save(k,11) = Linearfit(2);% intercept 
        k = k + 1;
    end
end
k = 1;
for i = max(Save(:,8)):-1:1
    for j = 1:length(Save(:,8))
        if Save(j,8) == i
            SavePoints1x(k,:) = Save(j - i + 1:j,1)';
            SavePoints1y(k,:) = Save(j - i + 1:j,2)';
            SavePoints2x(k,:) = Save(j - i + 1:j,3)';
            SavePoints2y(k,:) = Save(j - i + 1:j,4)';
            SaveTime(k,:) = Save(j - i + 1:j,5)';
            SavePos(k,:) = Save(j - i + 1:j,6)';
            SaveSlope(k,:) = Save(j - i + 1:j,7)';
            Save(j - i + 1:j,8) = 0;
            SaveMagPos(k,:) = Save(j - i + 1:j,9)';
            SaveSlopeFit(k,:) = Save(j - i + 1:j,10)';
            SaveInterceptFit(k,:) = Save(j - i + 1:j,11)';
            line([SavePoints1x(k,:) SavePoints2x(k,:)],[SavePoints1y(k,:) SavePoints2y(k,:)],'Linestyle','-','Linewidth',1,'color','g');hold on;
            line([SavePoints1x(k,:) SavePoints2x(k,:)],[SavePoints1x(k,:)*SaveSlopeFit(k,:) + SaveInterceptFit(k,:) SavePoints2x(k,:)*SaveSlopeFit(k,:) + SaveInterceptFit(k,:)],'Linestyle','-','Linewidth',1.5,'color','r');hold on;
            k = k + 1;
        end
    end
    k = 1;
    SAVE_Path = strcat('E:\Matlab analysis\',FileName_XYZ(1:end - length('.txt')),'_',num2str(Num_bead),'_bead_',num2str(i),'_Seg','.txt');
    fid = fopen(SAVE_Path,'wt');
    for l = 1:i
        fprintf(fid,'%s\t',strcat('point_1x_',num2str(l)));
    end
    for l = 1:i
        fprintf(fid,'%s\t',strcat('point_1y_',num2str(l)));
    end
    for l = 1:i
        fprintf(fid,'%s\t',strcat('point_2x_',num2str(l)));
    end
    for l = 1:i
        fprintf(fid,'%s\t',strcat('point_2y_',num2str(l)));
    end
    for l = 1:i
        fprintf(fid,'%s\t',strcat('time_',num2str(l)));
    end
    for l = 1:i
        fprintf(fid,'%s\t',strcat('position_',num2str(l)));
    end
    for l = 1:i
        fprintf(fid,'%s\t',strcat('slope_',num2str(l)));
    end
    for l = 1:i
        fprintf(fid,'%s\t',strcat('MagPos_',num2str(l)));
    end
    for l = 1:i
        fprintf(fid,'%s\t',strcat('slopeFit_',num2str(l)));
    end
    for l = 1:i
        fprintf(fid,'%s\t',strcat('InterceptFit_',num2str(l)));
    end
    fprintf(fid,'\n');
    matrix = [SavePoints1x SavePoints1y SavePoints2x SavePoints2y SaveTime SavePos SaveSlope SaveMagPos SaveSlopeFit SaveInterceptFit];
    [m_matrix,n_matrix] = size(matrix);
    for matrix_i = 1:m_matrix
        for matrix_j = 1:n_matrix
            if matrix_j == n_matrix
                fprintf(fid,'%g\n',matrix(matrix_i,matrix_j));
            else
                fprintf(fid,'%g\t',matrix(matrix_i,matrix_j));
            end
        end
    end
    fclose(fid);
    SavePoints1x = [];
    SavePoints1y = [];
    SavePoints2x = [];
    SavePoints2y = [];
    SaveTime = [];
    SavePos = [];
    SaveSlope = [];
    SaveMagPos = [];
    SaveSlopeFit = [];
    SaveInterceptFit = [];
end
delete Points.txt;
savefig(fig,strcat(SAVE_Path(1:end - 4),'.fig'));