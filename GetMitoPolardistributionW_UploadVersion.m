function [res,rr,grho0max]=GetMitoPolardistributionW_UploadVersion(im)
       
% analyze the polar distribution of mitochondrial image in cells
% input : im - gray scale image -
% ouput : res- polar intenity profile of mitochondrial image (N by 1
%        vector)
% `       rr - normalized radius vector (N by 1);
%        grho0max - value of max radius that used to normalize radius, rr
% 
% developed by Pei-hsun Wu @ Johns Hopkins University.

[yy,xx]=size(im);



[gx,gy]=meshgrid([1:xx],[1:yy]);

%         [bg,rbg]=estibkg(double(im),7);
bw0=im==0;
bw0=imopen(bw0,strel('disk',21));

bd=bwboundaries(~bw0);        
bw01=imclearborder(bw0); % inner circle map
grho0=bwdist(bw01); % get distance from inner circle;
st=regionprops(bw01,'Centroid'); % get centroid of inner circle
xyc=[st.Centroid];
xc=xyc(1); yc=xyc(2);
%         xyc=[xyc(1:2:end)' xyc(2:2:end)'];
[gtheta,~]=cart2pol(gx-xc,gy-yc); % obtain the orientation map
%         grho0=grho0*px/ra; % convert to um from px);
gthetad=gtheta*180/pi;
gthetad=round(gthetad);
resr=[];
for td=-180:1:180 % obtain the max r in a given orientation
    cc=gthetad==td;
    tmp2=grho0(cc & ~bw0);
    resr=[resr;[td,max(tmp2)]];
end
grho0max=interp1(resr(:,1),resr(:,2),gtheta/pi*180); % obtained max radius map at different angle
grho0n=grho0./grho0max; % obtained normalized radius map based on inner and outter circle                        

rr=[0:0.025:1.5]; % based on fixed size

[Irr,Irrs,N,Rrr,Irrsem]=deal(0);
for k=1:length(rr)-1
    cc=grho0n>rr(k) & grho0n<=rr(k+1);
    Irr(k)=mean(im(cc));
    Irrs(k)=std(single(im(cc)));
    N(k)=sum(cc(:));
    Rrr(k)=mean(grho0n(cc));
    Irrsem(k)=Irrs(k)./sqrt(N(k));
end
res(1:length(Irr),kp)=Irr(:);

%     ress=[ress;Irrs];
fprintf('\nfinish analyzing image');


return