close all;clear;
Num_bead = 41;
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
MagRos = [Motor(:,1) Motor(:,4)];
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

Smooth = 250;
Num_totalBead = fix((length(XYZ) - 5)/3);
Num_reference = Num_totalBead - 0;% 0 represents the last bead
framerate = 50;
for i = 5:3:5 + 3*Num_totalBead
    XYZsub{i} = - XYZ{i} + sgolayfilt(XYZ{5+ 3*Num_reference},4,Smooth - 1);
    %XYZsub{i} = - XYZ{i};
end
XYZSmooth = sgolayfilt(XYZsub{5 + Num_bead*3},4,Smooth - 1);

fig = figure('units','normalized','outerposition',[0 0.2 1 0.8]);
p = plot(XYZ{1}/framerate,XYZsub{5 + Num_bead*3},'.','MarkerSize',5,'hittest', 'off');hold on;
p = plot(XYZ{1}/framerate,XYZSmooth,'.','MarkerSize',5,'hittest', 'off');hold on;
set(gca,'position',[0.05 0.11 0.9 0.815]);
grid on;
grid minor;
ax = gca;
ax.GridAlpha = 0.7;
ax.MinorGridAlpha = 0.7;
title(strcat('DNA extension vs time-',FileName_XYZ(1:end - 4),'-bead-',num2str(Num_bead)));
xlabel('Time (s)');
ylabel('DNA extension (um)');
ylim([min(XYZSmooth(XYZSmooth ~= 0)) - 0.1,max(XYZSmooth(XYZSmooth ~= 0)) + 0.1]);
xlim([min(XYZ{1}/framerate) - 50, max(XYZ{1}/framerate) + 50]);
if exist('Points.txt', 'file')
    delete Points.txt;
end
set(gca,'buttondownfcn',@mybttnfcn_);

%open the Points file, calculate the duration and amplitude of each segment
%and save them.
[FileName_Save,PathName_Save] = uigetfile('D:\Program files\Matlab\myCode\.txt','Select a txt file');
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
        Save(k,5) = Save(k,3) - Save(k,1);
        Save(k,6) = mean(XYZ{5 + 3*Num_bead}(round(framerate*Points{1}(i)):round(framerate*Points{1}(i + 1))));% mean of segment extracted from raw data
        Save(k,7) = mean(XYZsub{5 + 3*Num_bead}(round(framerate*Points{1}(i)):round(framerate*Points{1}(i + 1))));% mean of segment extracted from raw data with dirft correction.
        Save(k,8) = Points{3}(i);
        Save(k,9) = MagRos(round(framerate*Points{1}(i)),2);% magnets position
        k = k + 1;
    end
end

% open hat curve file
[FileName_hat,PathName_hat] = uigetfile('D:\Program files\Matlab\myCode\.txt','Select a txt file');
[fid,message] = fopen(FileName_hat,'rt');
if fid == -1
    disp (message);
else
    disp (fid);
end
DELIMITER = '\t';
HEADERLINE = 1;
HatData = importdata(FileName_hat,DELIMITER,HEADERLINE);
fclose(fid);
HatCurve = HatData.data;
k = 1;
for i = max(Save(:,8)):-1:1
    for j = 1:length(Save(:,8))
        if Save(j,8) == i
            Save_points1x(k,:) = Save(j - i + 1:j,1)';
            Save_points1y(k,:) = Save(j - i + 1:j,2)';
            Save_points2x(k,:) = Save(j - i + 1:j,3)';
            Save_points2y(k,:) = Save(j - i + 1:j,4)';
            Save_time(k,:) = Save(j - i + 1:j,5)';
            Save_Pos(k,:) = Save(j - i + 1:j,6)';
            Save_PosDriftCorrected(k,:) = Save(j - i + 1:j,7)';
            Save(j - i + 1:j,8) = 0;
            Save_MagRos(k,:) = Save(j - i + 1:j,9)';
            line([Save_points1x(k,:) Save_points2x(k,:)],[Save_points1y(k,:) Save_points2y(k,:)],'Linestyle','-','Linewidth',1,'color','g');hold on;
            Save(j - i + 1:j,8) = 0;
            k = k + 1;
        end
    end
    if ~isempty(Save_MagRos)
        mMagRos = 1;
        mSlope = 1;
        for mMagRos = 1:length(Save_MagRos(:,1))
            for mHat = 1:length(HatCurve(:,Num_bead + 2))
                if round(Save_MagRos(mMagRos,1)) == round(HatCurve(mHat,1))
                    if round(Save_MagRos(mMagRos,1)) > 0
                        Save_HatSlope(mSlope,1) = HatCurve(mHat,Num_bead + 2) - HatCurve(mHat - 1,Num_bead + 2);
                        Save_HatSlope(mSlope,2) = HatCurve(mHat - 1,Num_bead + 2) - HatCurve(mHat - 2,Num_bead + 2);
                        Save_HatSlope(mSlope,3) = HatCurve(mHat - 2,Num_bead + 2) - HatCurve(mHat - 3,Num_bead + 2);
                        Save_HatSlope(mSlope,4) = mHat;
                        mSlope = mSlope + 1;
                    else
                        if round(Save_MagRos(mMagRos,1)) < 0
                            Save_HatSlope(mSlope,1) = HatCurve(mHat - 1,Num_bead + 2) - HatCurve(mHat,Num_bead + 2);
                            Save_HatSlope(mSlope,2) = HatCurve(mHat - 2,Num_bead + 2) - HatCurve(mHat - 1,Num_bead + 2);
                            Save_HatSlope(mSlope,3) = HatCurve(mHat - 3,Num_bead + 2) - HatCurve(mHat - 2,Num_bead + 2);
                            Save_HatSlope(mSlope,4) = mHat;
                            mSlope = mSlope + 1;
                        end
                    end
                end
            end
        end
    end
    if ~isempty(Save_MagRos)
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
            fprintf(fid,'%s\t',strcat('positionDriftCorrected_',num2str(l)));
        end
        fprintf(fid,'%s\t','MagRos');
        fprintf(fid,'%s\t','HatSlope1');
        fprintf(fid,'%s\t','HatSlope2');
        fprintf(fid,'%s\t','HatSlope3');
        fprintf(fid,'%s\t','HatSlopeTurn1');
        fprintf(fid,'\n');
        matrix = [Save_points1x Save_points1y Save_points2x Save_points2y Save_time Save_Pos Save_PosDriftCorrected Save_MagRos(:,1) Save_HatSlope];
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
        Save_points1x = [];
        Save_points1y = [];
        Save_points2x = [];
        Save_points2y = [];
        Save_time = [];
        Save_Pos = [];
        Save_PosDriftCorrected = [];
        Save_MagRos = [];
        Save_HatSlope = [];
    end
end
delete Points.txt
savefig(fig,strcat(SAVE_Path(1:end - 4),'.fig'));