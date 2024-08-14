function diaphragm_pos = THREE_diaphragm_motion(All_Ims,X,Y,force_use)

if nargin < 4
    force_use = 0;
end

%Function to evaluate all 2D images to get a measure of diaphragm position.
%Ultimately, it would be cool to include some of the neat 2D correlation
%matrix stuff, but for now just do a simple line (specified by X and Y)
%over the diaphragm.

%Start by making sure the line is in the right place
figure('Name','Check Line')
imagesc(All_Ims(:,:,1));
colormap(gray)
axis square
axis off
hold on
plot(X,Y,'g','LineWidth',2);
hold off

fail = 1;
while fail == 1
    Lines = zeros(max([X(2)-X(1),Y(2)-Y(1)]),size(All_Ims,3));
    dLines = zeros((max([X(2)-X(1),Y(2)-Y(1)]))-1,size(All_Ims,3));
    diaphragm_pos = zeros(1,size(All_Ims,3));
    for i = 1:size(All_Ims,3)
        pts = Lines(:,i);
        if Y(1) == Y(2)
            Lines(:,i) = All_Ims(Y(1):Y(2),(X(1)+1):X(2),i);
        else
            Lines(:,i) = All_Ims((Y(1)+1):Y(2),X(1):X(2),i);
        end
        diffline = abs(diff(Lines(:,i)));
        dLines(1:length(diffline),i) = abs(diff(Lines(:,i)));
        [~,diaphragm_pos(i)] = max(dLines(:,i));
    end
    %Let's look at Otsu thresholding as well.
    My_thresh = 0.5*graythresh(Lines);
    d_pos_2 = zeros(1,size(All_Ims,3));

    try
        if Lines(1,1)<Lines(end,1)
            for i = 1:size(All_Ims,3)
                d_pos_2(i) = find(Lines(:,i)<My_thresh,1,'last');
            end
        else
            for i = 1:size(All_Ims,3)
                d_pos_2(i) = find(Lines(:,i)>My_thresh,1,'last');
            end
        end
        fail = 0;
    catch
        disp('Otsu Mask Failed - Extending Line');
        if X(1) < X(2) 
            X(1) = X(1)-1;
            X(2) = X(2)+1;
        elseif X(1) > X(2)
            X(1) = X(1)+1;
            X(2) = X(2)-1;
        elseif Y(1) < Y(2) 
            Y(1) = Y(1)-1;
            Y(2) = Y(2)+1;
        elseif Y(1) > Y(2)
            Y(1) = X(1)+1;
            Y(2) = Y(2)-1;
        end
    end
    if force_use == 1
        fail = 0;
    end
end
figure('Name','Diaphragm Motion - No Line')
imagesc(Lines)
colormap(gray)
axis off

figure('Name','Diaphragm Motion - Otsu')
imagesc(Lines)
colormap(gray)
axis off
hold on
plot(1:size(All_Ims,3),d_pos_2,'r','LineWidth',1);
hold off
set(gcf,'Position',[37 475 1797 420]);
drawnow
figure('Name','Diaphragm Motion - Differential Method')
imagesc(Lines)
colormap(gray)
axis off
hold on
plot(1:size(All_Ims,3),diaphragm_pos,'g','LineWidth',1);
hold off
set(gcf,'Position',[37 475 1797 420]);

%% Sometimes Otsu's method works better - If so, uncomment this line
if force_use ~= 1 
    diaphragm_pos = d_pos_2;
end








