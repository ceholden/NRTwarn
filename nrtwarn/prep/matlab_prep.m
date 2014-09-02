MODISQA=hdfread('../../tests/MOD09G/MOD09GA.A2000055.h10v08.005.2006268014216.hdf','state_1km_1');

%% Mine Code
Water=( bitand(uint16(MODISQA),uint16(56))~=8 );
% Cloud Flag:  Cloud 0-2, Cirrus 8,9, Internal 10, Adjacent to cloud, 13

Cloud=mod(MODISQA,8) | mod(bitshift(MODISQA, -8),4) | mod(bitshift(MODISQA, -10),2) |  mod(bitshift(MODISQA, -13),2);
% Cloud=mod(MODISQA,8) | mod(bitshift(MODISQA, -8),4) | mod(bitshift(MODISQA, -10),2);
% Cloud=mod(MODISQA,8) | mod(bitshift(MODISQA, -10),2);
CloudBuffer=1-imerode(~Cloud,strel('square',7));
GoodLandObservation=(~CloudBuffer.*~Water);

%% Your Code
% temp=double(MODISQA);
% Mask_Land=floor(mod(temp+55,64)/63);
% Mask_Land(Mask_Land<1)=NaN;
% Mask_Cloud=mod(floor(temp/1024)+1,2).*floor(mod(temp+7,8)/7);
% Mask_Cloud(Mask_Cloud<1)=NaN;
% [I,J]=find(isnan(Mask_Cloud));
% if size(I,1)>0
% for Index=1:size(I,1)
% Mask_Cloud(max(I(Index)-3,1):min(I(Index)+3,size(Mask_Cloud,1)),max(J(Index)-3,1):min(J(Index)+3,size(Mask_Cloud,2)))=NaN;
% end
% end
% Mask=Mask_Land.*Mask_Cloud;

