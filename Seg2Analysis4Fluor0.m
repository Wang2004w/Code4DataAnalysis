close all;clear;
%open a .txt file
[FileName,PathName] = uigetfile('D:\Program files\Matlab\myCode\.txt','Select a txt file');
SAVEdata = [];
SAVEseq = [];
FileNameSave = [];
PointsOri1 = [];
PointsOri2 = [];
PointsOri3 = [];
Point1CH2 = [];
Point2CH2 = [];
Point3CH2 = [];
SegTime1 = [];
SegTime2 = [];
FRET1 = [];
FRET2 = [];
[fid,message] = fopen(FileName,'rt');
if fid == -1
    disp (message);
else
    disp (fid);
end
DELIMITER = '\t';
HEADERLINE = 1;
Fluorescence = importdata(FileName,DELIMITER,HEADERLINE);
fclose(fid);
if ~isempty(Fluorescence.data)
    Fluor1 = Fluorescence.data(1:2:end,:);
    Fluor2 = Fluorescence.data(2:2:end,:);
    if mean(Fluor1(1:10,2) - Fluor1(1:10,3)) > mean(Fluor2(1:10,2) - Fluor2(1:10,3))
        CH1 = Fluor1;
        CH2 = Fluor2;
    else
        if mean(Fluor1(1:10,2) - Fluor1(1:10,3)) < mean(Fluor2(1:10,2) - Fluor2(1:10,3))
            CH1 = Fluor2;
            CH2 = Fluor1;
        end
    end
    xCH1 = CH1(:,1);
    xCH2 = CH2(:,1);
    IDD = CH1(:,2) - CH1(:,3);
    IDA = CH1(:,4) - CH1(:,5);
    IAA = CH2(:,4) - CH2(:,5);
    fig = figure('units','normalized','outerposition',[0 0.4 1 0.6]);
    pIDD = plot(xCH1,IDD,'.-','LineWidth',1,'MarkerSize',12,'color',[0.4660 0.6740 0.1880],'hittest', 'off');hold on;
    %pIDA = plot(xCH1,IDA,'.-','LineWidth',1,'MarkerSize',12,'color',[0.8500 0.3250 0.0980],'hittest', 'off');hold on;
    pIAA = plot(xCH2,IAA,'.-','LineWidth',1,'MarkerSize',12,'color',[0.4940 0.1840 0.5560],'hittest', 'off');hold on;
    line(xlim(),[0,0],'LineWidth',1.5,'color','k');
    grid on;
    grid minor;
    ax = gca;
    ax.GridAlpha = 0.7;
    ax.MinorGridAlpha = 0.7;
    title(FileName(1:end - 4));
    xlabel('Time (ms)');
    ylabel('Intensity (a.u.)');
    IDDDAAA = [IDD;IDA;IAA];
    ylim([min(IDDDAAA(IDDDAAA ~= 0)) - 0.1,max(IDDDAAA(IDDDAAA ~= 0)) + 0.1]);
    xlim([min(xCH1) - 10, max(xCH1 + 10)]);
    if exist('Points.txt', 'file')
        delete Points.txt;
    end
    set(gca,'buttondownfcn',@mybttnfcn_);
    [FileNamePoints,PathNamePoints] = uigetfile('D:\Program files\Matlab\myCode\MagneticTrapDataAnalysis\.txt','Select a txt file');
    [fid,message] = fopen(FileNamePoints,'rt');
    if fid == -1
        disp (message);
    else
        disp (fid);
    end
    Points = textscan(fid,'','delimiter','\t','HeaderLines',0,'TreatAsEmpty',{'1.#QNAN00'});
    fclose(fid);
    for i = 2:length(xCH2) - 1
        if abs(Points{1}(1) - xCH2(i)) <= abs(Points{1}(1) - xCH2(i - 1))&&abs(Points{1}(1) - xCH2(i)) <= abs(Points{1}(1) - xCH2(i + 1))
            Point1CH2(1,1) = xCH2(i);
            Point1SeqCH2 = i;
        end
        if abs(Points{1}(2) - xCH2(i)) <= abs(Points{1}(2) - xCH2(i - 1))&&abs(Points{1}(2) - xCH2(i)) <= abs(Points{1}(2) - xCH2(i + 1))
            Point2CH2(1,1) = xCH2(i);
            Point2SeqCH2 = i;
        end
        if abs(Points{1}(3) - xCH2(i)) <= abs(Points{1}(3) - xCH2(i - 1))&&abs(Points{1}(3) - xCH2(i)) <= abs(Points{1}(3) - xCH2(i + 1))
            Point3CH2(1,1) = xCH2(i);
            Point3SeqCH2 = i;
        end
    end
    FileNameSave(1,1) = str2double(FileName(1:end - length('.txt')));
    PointsOri1(1,1) = Points{1}(1);
    PointsOri2(1,1) = Points{1}(2);
    PointsOri3(1,1) = Points{1}(3);
    SegTime1(1,1) = Point2CH2(1,1) - Point1CH2(1,1);
    SegTime2(1,1) = Point3CH2(1,1) - Point2CH2(1,1);
    FRET1(1,1) = mean(IDA(Point1SeqCH2:Point2SeqCH2))/(mean(IDA(Point1SeqCH2:Point2SeqCH2)) + mean(IDD(Point1SeqCH2:Point2SeqCH2)));
    FRET2(1,1) = mean(IDA(Point2SeqCH2:Point3SeqCH2))/(mean(IDA(Point2SeqCH2:Point3SeqCH2)) + mean(IDD(Point2SeqCH2:Point3SeqCH2)));
    if exist('Points.txt', 'file')
        delete Points.txt;
    end
    SAVEPath = strcat('E:\Matlab analysis\',FileName(1:end - length('.txt')),'_SegAnalysis4Fluor');
    savefig(fig,strcat(SAVEPath,'.fig'));
    close Figure 1;
end
SAVEPath = strcat('E:\Matlab analysis\',FileName(1:end - length('.txt')),'_SegAnalysis4Fluor.txt');
fid = fopen(SAVEPath,'wt');
fprintf(fid,'%s\t','FileNameSeq');
fprintf(fid,'%s\t','PointsOri1');
fprintf(fid,'%s\t','PointsOri2');
fprintf(fid,'%s\t','PointsOri3');
fprintf(fid,'%s\t','Point1CH2');
fprintf(fid,'%s\t','Point2CH2');
fprintf(fid,'%s\t','Point3CH2');
fprintf(fid,'%s\t','SegTime1');
fprintf(fid,'%s\t','SegTime2');
fprintf(fid,'%s\t','FRET1');
fprintf(fid,'%s\t','FRET2');
fprintf(fid,'\n');
matrix = [FileNameSave PointsOri1 PointsOri2 PointsOri3 Point1CH2 Point2CH2 Point3CH2 SegTime1 SegTime2 FRET1 FRET2];
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