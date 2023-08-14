function mybttnfcn_(h,~)
hf = get(h,'parent');
b = get(hf,'selectiontype');
cp = get(gca,'CurrentPoint');
xy = getappdata(hf,'xypoints');
xy_leftClick = getappdata(hf,'xypoints_leftClick');
if strcmpi(b,'normal')
    text(cp(1,1),cp(1,2),'')
    line(cp(1,1),cp(1,2),'linestyle','none','marker','o','LineWidth',1.5,'color','k');
    xy = [xy;cp(1,1:2)];
    xy_leftClick = [xy_leftClick;cp(1,1:2)];
    setappdata(hf,'xypoints',xy);
    setappdata(hf,'xypoints_leftClick',xy_leftClick);
elseif strcmpi(b,'alt')
    %text(cp(1,1),cp(1,2),'Saved')
    xy = getappdata(hf,'xypoints');
    xy_leftClick = getappdata(hf,'xypoints_leftClick');
    XY = [xy_leftClick (1:length(xy_leftClick(:,1)))'];
    line(xy_leftClick(:,1),xy_leftClick(:,2),'LineWidth',1,'color','k')
    xy_leftClick = [];
    setappdata(hf,'xypoints_leftClick',xy_leftClick);
    fid = fopen('D:\Program Files\MATLAB\R2017a\bin\myCode\FluorescenceAnalysis\Points.txt','at');
    matrix = XY;%´ýÊä³ö¾ØÕó
    [m_matrix,n_matrix] = size(matrix);
    for i = 1:1:m_matrix
        for j = 1:1:n_matrix
            if j == n_matrix
                fprintf(fid,'%g\n',matrix(i,j));
            else
                fprintf(fid,'%g\t',matrix(i,j));
            end
        end
    end
    fclose(fid);
else
    text(xy(1,1),xy(1,2),'Careful there, crazy man!')    
end
