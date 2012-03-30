clc
clear all
load OutputLog.dat
% Loop through each trajectory, changing the data points.

n=size(OutputLog,2);
t=[0:0.1:2*pi];
radius=0.1;
dt=0.05;
test=0;
for k=1:size(OutputLog,1)-1
    for i=2:3:n
        hold on
        axis([-1 1 -1 1])
	title(num2str(OutputLog(k,1))) 
        switch test
            case 2;
                %line([OutputLog(k,i) OutputLog(k+1,i)],[OutputLog(k,i+1) OutputLog(k+1,i+1)],'color','b'); 
                x=radius*sin(t)+OutputLog(k,i);
                y=radius*cos(t)+OutputLog(k,i+1);
                plot(x,y,'color','b')
            case 4;
                x=radius*sin(t)+OutputLog(k,i);
                y=radius*cos(t)+OutputLog(k,i+1);
                plot(x,y,'color','r')
                %line([OutputLog(k,i) OutputLog(k+1,i)],[OutputLog(k,i+1) OutputLog(k+1,i+1)],'color','r');   
            case 6;
                x=radius*sin(t)+OutputLog(k,i);
                y=radius*cos(t)+OutputLog(k,i+1);
                plot(x,y,'color','g')
                %line([OutputLog(k,i) OutputLog(k+1,i)],[OutputLog(k,i+1) OutputLog(k+1,i+1)],'color','g');  
            case 8;
                x=radius*sin(t)+OutputLog(k,i);
                y=radius*cos(t)+OutputLog(k,i+1);
                plot(x,y,'color','k')
                %line([OutputLog(k,i) OutputLog(k+1,i)],[OutputLog(k,i+1) OutputLog(k+1,i+1)],'color','k'); 
            case 10;
                x=radius*sin(t)+OutputLog(k,i);
                y=radius*cos(t)+OutputLog(k,i+1);
                plot(x,y,'color','m') 
           case 12;
                x=radius*sin(t)+OutputLog(k,i);
                y=radius*cos(t)+OutputLog(k,i+1);
                plot(x,y,'color','c')
	   otherwise;
                x=radius*sin(t)+OutputLog(k,i);
                y=radius*cos(t)+OutputLog(k,i+1);
                plot(x,y,'color','k')

        end
    end
    pause(dt);
    clf
end




