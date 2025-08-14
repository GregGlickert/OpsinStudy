close all;
clear all;
load('TransmittanceData');%%%%Drew's transmittancedata
figure(1)
surf(x,y,transImg,'EdgeColor','interp','FaceColor','interp','FaceLighting','phong')
title('transmittance');colorbar
location=load('location.txt')*1e3 %%%load cell location,um
upscale=27;
PCELL=900*upscale;
INTCELL=100*upscale;   %%according to Drew, only 10% is FS.
NCELL=PCELL+INTCELL;
figure(2)
x_p=location(1:PCELL,1);y_p=location(1:PCELL,2);z_p=location(1:PCELL,3);
x_I=location(PCELL+1:NCELL,1);y_I=location(PCELL+1:NCELL,2);z_I=location(PCELL+1:NCELL,3);
scatter3(x_p,y_p,z_p,'blue','filled');hold on;
scatter3(x_I,y_I,z_I,'filled','red');hold on;

%%%%imagine light is shed from the top center of structure