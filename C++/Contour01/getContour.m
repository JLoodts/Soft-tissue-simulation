clear all;
close all;
% fotootje laden
RGB = imread('onion.jpg','jpg');
X = rgb2gray(RGB);
imagesc(X);
colormap gray;
axis image;


axis manual
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
		if( length(x) > 1 )
			plot(x([end-1 end]),y([end-1 end]),'r-');
		end
		plot(px,py,'ro');
	end
end

hold off;
 
% load clickResult.mat;

coo = [x;y];
min = min(coo,[],2);
coo = coo-min*ones(1,length(x));
range = max(max(coo,[],2),[],1);
coo = coo./(range*ones(2,length(x)))
x = coo(1,:);
y = -coo(2,:)+ones(1,length(y));

save clickResult.mat x y;
%save cornerPoints.txt coo -ASCII -DOUBLE

