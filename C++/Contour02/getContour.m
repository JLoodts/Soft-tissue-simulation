rangex = [0 9.23e-3]; % real height of the entire picture
rangey = [0 1.77e-3]; % real width of the entire picture
dx = 1.5*rangey(2); % part of the image in the x-direction that is visible

RGB = imread('onion.jpg','jpg');
imshow(rangex,rangey,RGB);
iptsetpref('ImshowAxesVisible', 'on');
iptsetpref('ImshowBorder','loose');
axis([rangex, rangey]);
a = gca;
%%%%% Set appropriate axis limits and settings
set(gcf,'doublebuffer','on');
% This avoids flickering when updating the axis
set(a,'xlim',[0 1e-3]);
set(a,'ylim',[min(rangey) max(rangey)]);

%%%%% Generate constants for use in uicontrol initialization
pos=get(a,'position');
Newpos=[pos(1) pos(2)-0.1 pos(3) 0.05];
% This will create a slider which is just underneath the axis
% but still leaves room for the axis labels above the slider
xmax=max(rangex);
S=['set(gca,''xlim'',get(gcbo,''value'')+[0 ' num2str(dx) '])'];
% Setting up callback string to modify XLim of axis (gca)
% based on the position of the slider (gcbo)

%%%%% Creating Uicontrol
h=uicontrol('style','slider',...
    'units','normalized','position',Newpos,...
    'callback',S,'min',0,'max',xmax-dx);

title('Click left to draw polyline, click right to terminate')
hold on;

% herhaal tot andere dan linkermuisknop ingedrukt
x = []; y = [];
while(1)
	[px,py,button] = ginput(1);
	if( button ~= 1 )
		break;
	else
		x = [x px] ; y = [y py];
% 		if( length(x) > 1 )
% 			plot(x([end-1 end]),y([end-1 end]),'r-');
% 		end
		plot(px,py,'ro');   
    end
end

hold off;
 
% load clickResult.mat;

% coo = [x;y];
% min = min(coo,[],2);
% coo = coo-min*ones(1,length(x));
% range = max(max(coo,[],2),[],1);
% coo = coo./(range*ones(2,length(x)))
% x = coo(1,:);
% y = -coo(2,:)+ones(1,length(y));
% % herschalen naar de grootste dimensie (in dit geval de lengte van de strip)
% x = x*8.91e-3;
% y = y*8.91e-3;
ymax = max(rangey);
y = ymax - y;

save clickResult.mat x y;
fid = fopen('cornerPoints.txt','w');
fprintf(fid,'%3.0f %12.6f  %12.6f\n',[[1:1:length(x)]; x; y]);
fclose(fid);
%save cornerPoints.txt coo -ASCII -DOUBLE