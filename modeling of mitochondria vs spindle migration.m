% Program for “Dynamic Organelle Distribution Initiates Actin-Based 
%              Spindle Migration in Mouse Oocytes”
% Contact: Yizeng Li (liyizeng52@hotmail.com)
%          Sean X. Sun (ssun@jhu.edu)

clear
clc

Umt0 = 0.12d-18;
Umtp = 0.5*Umt0;
dmt0p = 0.13d-6;

Ums0 = 1.25d-18;
Umsp = 1.2d-18;
dms0p = 0.1d-6;


as = 15.d-6;     % (m) long axis of the spindle
bs = 11.d-6;      % (m) short axis of the spindle
Rmt = 0.15d-6;    % (m) radius of the Mitos
rmt = 0.45d-6;    % (m) initial distance between the centers of Mitos to the boundary of the spindle

width = 5.d-6;      % (m) half bandwidth to calculate the number of Mito

Layer = 8;
Nfa = 160*2;        % number of actin fialment around the spindle
nmt = 180;
Nmt = nmt*Layer;       % number of Mitos
dt = 1;             % (s) time step
Nt = 210*60/dt+1;        % number of time steps
Time = dt:dt:dt*Nt; % (s) time vector


etasa = 4d-4;%2*13.0d-3;    % (N s/m) drag coefficient along the long axis of the spindle
etasb = 8d-4;%2*35.0d-3;    % (N s/m) drag coefficient along the short axis of the spindle
etasr = 0.1d-4;    % (N s m) drag coefficient of rotation of the spindle
etamt = 3d-6;%1.5d-3;     % (N s/m) drag coefficient of Miots. Control: 1.5d-3; increased: 2.5d-3

LfaMax = 0.35d-6;     % (m) maximum lenght of an actin filament
FfaMax = 1.d-11;    % (N) maximum force in an actin filament. Control: 1.d-11; reduced: 1.d-12;
damc = 0.05d-6;    % (m) critical distance of potential forces between Mitos and actin filaments
kma = 0.5d-4;         % coefficient of potential force between filaments and Mito
Fon = 50.d-12;    % (N) cut-off force for polymerization
Foff = 50.d-12;   % (N) cut-off force for deploymerization
Kon = 10;       % max. rate of polymerization
Koff = 3;       % max. rate of depolymerization
delta = 5.d-9;  % (m) size of a G-actin

dmt0 = dmt0p*(2+sqrt(4-3*Umtp/Umt0))/3;
bmt = 2*Umt0*((dmt0/dmt0p)^2-dmt0/dmt0p)*dmt0;
Frep0 = 2*Umt0/dmt0p;
Fatr0 = -(-bmt/dmt0^2);

sigma = 6d-16;    % (N) band with of thermal noise

dmt0s = dms0p*(2+sqrt(4-3*Umsp/Ums0))/3;
bms = 2*Ums0*((dmt0s/dms0p)^2-dmt0s/dms0p)*dmt0s;
Frep0s = 2*Ums0/dms0p;
Fatr0s = -(-bms/dmt0s^2);

Xs = zeros(Nt,2);
Thetas = zeros(Nt,1);

Xsa = zeros(Nt, Nfa);  % point of the actin on the spindle
Ysa = zeros(Nt, Nfa);
Xam = zeros(Nt, Nfa); % the other end of the actin filaments close to Mitos
Yam = zeros(Nt, Nfa);
Xmt = NaN(Nt, Nmt);
Ymt = NaN(Nt, Nmt);
Xsm = zeros(Nt, Nmt); % closest point on the spindle to each Mito
Ysm = zeros(Nt, Nmt); 

Lfa = zeros(Nt, Nfa);
Dam = NaN(Nt, Nfa);  % distance (directional) between the tip of the actin filament to the Mitos
Fam = zeros(Nt, Nfa); % potential force between the tip of actin filaments and the Mitos
Link = zeros(Nt, Nfa); % matric that store information on the force relation between actin and Mito

SMR = zeros(Nt,1);

%% Initial Distribution
tic

thetaDist = rand(1,Nfa)*2*pi;
RDist = as*bs./((bs*cos(thetaDist)).^2+(as*sin(thetaDist)).^2).^(1/2);
Xsa(1,:) = Xs(1,1) + RDist.*cos(thetaDist+Thetas(1));
Ysa(1,:) = Xs(1,2) + RDist.*sin(thetaDist+Thetas(1));

Lfa(1,:) = LfaMax*(0.5+0.5*rand(1,Nfa));
thetafa = thetaDist + (2*pi/Nfa)*(rand(1,Nfa)-0.5);
Xam(1,:) = Xsa(1,:) + Lfa(1,:).*cos(thetafa+Thetas(1));
Yam(1,:) = Ysa(1,:) + Lfa(1,:).*sin(thetafa+Thetas(1));

thetaMt = linspace(0,2*pi,nmt+1);
thetaMt = thetaMt(1:nmt);
RDist = as*bs./((bs*cos(thetaMt)).^2+(as*sin(thetaMt)).^2).^(1/2);
dtheta = 2*pi/nmt/2;
for iL = 1:Layer
    Xmt(1,(iL-1)*round(Nmt/Layer)+1:iL*nmt) = Xs(1,1)+(RDist+rmt+(iL-1)*2.5*Rmt).*cos(thetaMt+Thetas(1)+dtheta*mod(iL,2));
    Ymt(1,(iL-1)*round(Nmt/Layer)+1:iL*nmt) = Xs(1,2)+(RDist+rmt+(iL-1)*2.5*Rmt).*sin(thetaMt+Thetas(1)+dtheta*mod(iL,2));
end

thetaplot = linspace(0,2*pi,401);
rplot = as*bs./((bs*cos(thetaplot)).^2+(as*sin(thetaplot)).^2).^(1/2);

figure(1)
plot((Xs(1,1)+rplot.*cos(thetaplot+Thetas(1)))*1d6, (Xs(1,2)+rplot.*sin(thetaplot+Thetas(1)))*1d6,'-r','linewidth',3); hold on
for i1 = 1:Nfa
    plot([Xsa(1,i1),Xam(1,i1)]*1d6,[Ysa(1,i1),Yam(1,i1)]*1d6,'-k','linewidth',2);
end
for i1 = 1:Nmt
    plot((Xmt(1,i1)+Rmt*cos(thetaplot))*1d6, (Ymt(1,i1)+Rmt*sin(thetaplot))*1d6,'-','Color',[0,0.5,0],'linewidth',2);
end
hold off
set(gca, 'fontsize',13);
axis equal
axis([-35 35 -25 25]);
set(gcf,'Color','w');
axis off

%% Time Step

Indframe = 0;
for i = 1:Nt-1
    
   
    %% Get the distance and force between the tip of each Mito and the directional closest Mito
    for i1 = 1:Nfa
        replace = 1;
        while replace == 1
           temp = ((Xam(i,i1)-Xmt(i,:)).^2+(Yam(i,i1)-Ymt(i,:)).^2).^(1/2);
           if min(temp) > Rmt && Lfa(i,i1) < LfaMax
               replace = 0;
           else
               Lfa(i,i1) = LfaMax*(0.5*0+2*0.5*rand(1));
               
               rand_temp = rand(1);
               if rand_temp < 0.5
                   thetaDist(i1) = rand(1)*2*pi;
               else
                   if i < Nt/1
                       thetaDist(i1) = normrnd(pi,0.4,1,1);
                   else
                       thetaDist(i1) = normrnd(0,0.4,1,1);
                   end
               end
                   
               thetafa(i1) = thetaDist(i1) + (2*pi/Nfa)*(rand(1)-0.5);
               RDist(i1) = as*bs./((bs*cos(thetaDist(i1))).^2+(as*sin(thetaDist(i1))).^2).^(1/2);
               Xsa(i,i1) = Xs(i,1) + RDist(i1).*cos(thetaDist(i1)+Thetas(i));
               Ysa(i,i1) = Xs(i,2) + RDist(i1).*sin(thetaDist(i1)+Thetas(i));
               Xam(i,i1) = Xsa(i,i1) + Lfa(i,i1).*cos(thetafa(i1)+Thetas(i));
               Yam(i,i1) = Ysa(i,i1) + Lfa(i,i1).*sin(thetafa(i1)+Thetas(i));
           end
        end
        Lfa(i,i1) = ((Xsa(i,i1)-Xam(i,i1)).^2 + (Ysa(i,i1)-Yam(i,i1)).^2).^(1/2);
        nfa = [(Xam(i,i1)-Xsa(i,i1))./Lfa(i,i1); (Yam(i,i1)-Ysa(i,i1))./Lfa(i,i1)];
        
        XmtTemp = Xmt(i,:);
        YmtTemp = Ymt(i,:);
        
        Find = 0;
        ifind = 1;
        while Find == 0 && ifind < 3 % Nmt
            DTemp = ((Xam(i,i1)-XmtTemp).^2+(Yam(i,i1)-YmtTemp).^2).^(1/2);
            [MinTemp, iTemp] = min(DTemp);
            Dsamt1 =  ((XmtTemp(iTemp)-Xsa(i,i1)).^2 + (YmtTemp(iTemp)-Ysa(i,i1)).^2).^(1/2);
            nsamt = [(XmtTemp(iTemp)-Xsa(i,i1))/Dsamt1; (YmtTemp(iTemp)-Ysa(i,i1))/Dsamt1];
            cosTemp = nfa(1)*nsamt(1) + nfa(2)*nsamt(2);
            sinTemp = (1-cosTemp^2)^(1/2);
            Dsamt2 = Dsamt1*sinTemp;
            if Dsamt2 > Rmt || cosTemp < 0
                XmtTemp(iTemp) = NaN;
                YmtTemp(iTemp) = NaN;
            else
                Daemt3 = (Rmt^2-Dsamt2^2)^(1/2);
                Dam(i,i1) = Dsamt1*cosTemp - Lfa(i,i1) - Daemt3;
                Link(i,i1) = iTemp;
                Find = 1;
            end
            ifind = ifind + 1;
        end
        if Dam(i,i1) < damc
            Fam(i,i1) = kma*(damc - Dam(i,i1));
        else
            Fam(i,i1) = 0;
        end
    end
    
    %% Get the distance and force between each pair of Mitos and Mitos to the spindle
    Dmt = zeros(Nmt+1,Nmt+1);
    Fmt = zeros(Nmt+1,Nmt+1);
    Xspindle = Xs(i,1)+rplot.*cos(thetaplot+Thetas(i));
    Yspindle = Xs(i,2)+rplot.*sin(thetaplot+Thetas(i));
    Thermal = normrnd(0,sigma,Nmt-1,1);
    for imt = 1:Nmt-1
        % With other Mitos
        Dmt(imt,imt+1:Nmt) = ((Xmt(i,imt)-Xmt(i,imt+1:Nmt)).^2+(Ymt(i,imt)-Ymt(i,imt+1:Nmt)).^2).^(1/2)-2*Rmt;
        for imt2 = imt+1:Nmt
            if Dmt(imt,imt2) > dmt0
                Fmt(imt,imt2) = -Fatr0*dmt0^2/Dmt(imt,imt2)^2 + Thermal(imt);
            else
                Fmt(imt,imt2) = -(Frep0+Fatr0)/dmt0*Dmt(imt,imt2) + Frep0 + Thermal(imt);
            end
        end
        % With the spindle
        DAll = ((Xmt(i,imt)-Xspindle).^2+(Ymt(i,imt)-Yspindle).^2).^(1/2)-Rmt;
        [Dmt(imt,Nmt+1), iS] = min(DAll);
        Xsm(i,imt) = Xspindle(iS);
        Ysm(i,imt) = Yspindle(iS);
        if Dmt(imt,Nmt+1) > dmt0s
            Fmt(imt,Nmt+1) = -Fatr0s*dmt0s^2/Dmt(imt,Nmt+1)^2;
        else
            Fmt(imt,Nmt+1) = -(Frep0s+Fatr0s)/dmt0s*Dmt(imt,Nmt+1) + Frep0s;
        end
    end
    % With the spindle (for the last Mito)
    for imt = Nmt:Nmt
        DAll = ((Xmt(i,imt)-Xspindle).^2+(Ymt(i,imt)-Yspindle).^2).^(1/2)-Rmt;
        [Dmt(imt,Nmt+1), iS] = min(DAll);
        Xsm(i,imt) = Xspindle(iS);
        Ysm(i,imt) = Yspindle(iS);
        if Dmt(imt,Nmt+1) > dmt0s
            Fmt(imt,Nmt+1) = -Fatr0s*dmt0s^2/Dmt(imt,Nmt+1)^2;
        elseif Dmt(imt,Nmt+1) < dmt0s
            Fmt(imt,Nmt+1) = -(Frep0s+Fatr0s)/dmt0s*Dmt(imt,Nmt+1) + Frep0s;
        end
    end
    Dmt = Dmt + Dmt';
    Fmt = Fmt + Fmt';
    
    %% plot
    if mod(i,round(120/dt)) == 0
        fprintf('Time Step: %d\n',i);
        Indframe = Indframe + 1;
        toc
        
        figure(2)
        plot((Xs(i,1)+rplot.*cos(thetaplot+Thetas(i)))*1d6, (Xs(i,2)+rplot.*sin(thetaplot+Thetas(i)))*1d6,'-r','linewidth',3); hold on
        for i1 = 1:Nfa
            plot([Xsa(i,i1),Xam(i,i1)]*1d6,[Ysa(i,i1),Yam(i,i1)]*1d6,'-k','linewidth',2);
        end
        for i1 = 1:Nmt
            rectangle('Position',[Xmt(i,i1)-Rmt Ymt(i,i1)-Rmt 2*Rmt 2*Rmt]*1d6,...
                'Curvature',1,'FaceColor',[0,0.5,0],'EdgeColor',[0,0.5,0],'LineWidth',0.1);
        end
        hold off
        set(gca, 'fontsize',13);
        xlabel('x ({\mu}m)', 'fontsize',13);
        ylabel('y ({\mu}m)', 'fontsize',13);
        title([num2str(i*dt/60),' min'],'FontSize',13);
        axis equal
        axis([-45 45 -35 35])
        set(gcf,'Color','w');
        mov(Indframe) = getframe(gcf);
        pause(0.01)
        
    end
    %% calculation of the motion of the spindle
    
    norma = [cos(Thetas(i)),sin(Thetas(i))];
    normb = [-sin(Thetas(i)),cos(Thetas(i))];
    unit_fa = [cos(thetafa+Thetas(i))',sin(thetafa+Thetas(i))'];
    unit_mt = -[(Xmt(i,:)-Xsm(i,:))'./(Dmt(1:Nmt,Nmt+1)+Rmt), (Ymt(i,:)-Ysm(i,:))'./(Dmt(1:Nmt,Nmt+1)+Rmt)];
    
    fa = 0;
    fb = 0;
    fr = 0;
    
    fa = fa - sum(Fam(i,:)'.*(unit_fa(:,1)*norma(1) + unit_fa(:,2).*norma(2)));
    fb = fb - sum(Fam(i,:)'.*(unit_fa(:,1)*normb(1) + unit_fa(:,2).*normb(2)));
    rs = [Xsa(i,:)'-Xs(i,1), Ysa(i,:)'-Xs(i,2)];
    fr = fr - sum(Fam(i,:)'.*(rs(:,1).*unit_fa(:,2)-rs(:,2).*unit_fa(:,1)));
       
    fa = fa + sum(Fmt(1:Nmt,Nmt+1).*(unit_mt(:,1)*norma(1) + unit_mt(:,2)*norma(2)));
    fb = fb + sum(Fmt(1:Nmt,Nmt+1).*(unit_mt(:,1)*normb(1) + unit_mt(:,2)*normb(2)));
    rs = [Xsm(i,:)'-Xs(i,1), Ysm(i,:)'-Xs(i,2)];
    fr = fr + sum(Fmt(1:Nmt,Nmt+1).*(rs(:,1).*unit_mt(:,2)-rs(:,2).*unit_mt(:,1)));
    
    dXs_a = fa/etasa;
    dXs_b = fb/etasb;
    dXs = dXs_a*norma(1) + dXs_b*normb(1);
    dYs = dXs_a*norma(2) + dXs_b*normb(2);
    dXs_abs = (dXs^2+dYs^2)^(1/2);
    dThetas = fr/etasr;
    
    Xs(i+1,1) = Xs(i,1) + dt*dXs;
    Xs(i+1,2) = Xs(i,2) + dt*dYs;
    Thetas(i+1) = Thetas(i) + dt*dThetas;
    
    %% calcuation of the motion of the actin filaments
    kon  = -Kon/Fon*Fam(i,:) + Kon;
    koff = Koff/Foff*Fam(i,:);
    dLfa = (kon-koff)*delta;
    Lfa(i+1,:) = Lfa(i,:) + dLfa;
    
    RDist = as*bs./((bs*cos(thetaDist)).^2+(as*sin(thetaDist)).^2).^(1/2);
    Xsa(i+1,:) = Xs(i+1,1) + RDist.*cos(thetaDist+Thetas(i+1));
    Ysa(i+1,:) = Xs(i+1,2) + RDist.*sin(thetaDist+Thetas(i+1));
    Xam(i+1,:) = Xsa(i+1,:) + Lfa(i+1,:).*cos(thetafa+Thetas(i+1));
    Yam(i+1,:) = Ysa(i+1,:) + Lfa(i+1,:).*sin(thetafa+Thetas(i+1));
    
    %% calcuation of the motion of the Mitos
    for imt = 1:Nmt
        fmtx = 0;
        fmty = 0;
        mask = Link(i,:) == imt;
        indexf = find(mask);
        fmtx = fmtx + sum(unit_fa(indexf,1)'.*Fam(i,indexf));
        fmty = fmty + sum(unit_fa(indexf,2)'.*Fam(i,indexf));
        
        unit = [(Xmt(i,imt)-Xmt(i,1:Nmt))'./(Dmt(imt,1:Nmt)+2*Rmt)', (Ymt(i,imt)-Ymt(i,1:Nmt))'./(Dmt(imt,1:Nmt)+2*Rmt)'];
        fmtx = fmtx + sum(unit(:,1)'.*Fmt(imt,1:Nmt));
        fmty = fmty + sum(unit(:,2)'.*Fmt(imt,1:Nmt));
            
        unit = [(Xmt(i,imt)-Xsm(i,imt)), (Ymt(i,imt)-Ysm(i,imt))]/(Dmt(imt,Nmt+1)+Rmt); 
        fmtx = fmtx + unit(1)*Fmt(imt,Nmt+1);
        fmty = fmty + unit(2)*Fmt(imt,Nmt+1);
        
        dmtx = fmtx/etamt;
        dmty = fmty/etamt;
        Xmt(i+1,imt) = Xmt(i,imt) + dt*dmtx;
        Ymt(i+1,imt) = Ymt(i,imt) + dt*dmty;
    end
        
    
    %% Check if any Mito goes into the spindle
    dtemp = ((Xmt(i+1,:)-Xs(i+1,1)).*cos(Thetas(i+1))+(Ymt(i+1,:)-Xs(i+1,2)).*sin(Thetas(i+1))).^2/as^2 ...
        + ((Xmt(i+1,:)-Xs(i+1,1)).*sin(Thetas(i+1))-(Ymt(i+1,:)-Xs(i+1,2)).*cos(Thetas(i+1))).^2/bs^2;
    mask = dtemp < 1;
    
    vect = [Xmt(i+1,:) - Xs(i+1,1); Ymt(i+1,:) - Xs(i+1,2)];
    dist = (vect(1,:).^2 + vect(2,:).^2).^(1/2);
    norm = [vect(1,:)./dist; vect(2,:)./dist];
    theta = get_angle(norm(1,:),norm(2,:));
    RDist = as*bs./((bs*cos(theta)).^2+(as*sin(theta)).^2).^(1/2);
    
    Xmt(i+1,mask) = Xs(i+1,1) + (RDist(mask)+2*Rmt).*cos(theta(mask));
    Ymt(i+1,mask) = Xs(i+1,2) + (RDist(mask)+2*Rmt).*sin(theta(mask));
       
    SMR(i+1) = SMR(i) + (dXs_abs^2)*dt^2;
end

%% Reduce data size

SaveTime = 60/dt;  % sec

Xs = Xs(1:SaveTime:Nt,:);
Thetas = Thetas(1:SaveTime:Nt);

Xsa = Xsa(1:SaveTime:Nt, :);  % point of the actin on the spindle
Ysa = Ysa(1:SaveTime:Nt, :);
Xam = Xam(1:SaveTime:Nt, :); % the other end of the actin filaments close to Mitos
Yam = Yam(1:SaveTime:Nt, :);
Xmt = Xmt(1:SaveTime:Nt, :);
Ymt = Ymt(1:SaveTime:Nt, :);
Xsm = Xsm(1:SaveTime:Nt, :); % closest point on the spindle to each Mito
Ysm = Ysm(1:SaveTime:Nt, :); 

Lfa = Lfa(1:SaveTime:Nt, :);
Dam = Dam(1:SaveTime:Nt, :);  % distance (directional) between the tip of the actin filament to the Mitos
Fam = Fam(1:SaveTime:Nt, :); % potential force between the tip of actin filaments and the Mitos
Link = Link(1:SaveTime:Nt, :); % matric that store information on the force relation between actin and Mito

SMR = SMR(1:SaveTime:Nt);

Nt2 = length(Thetas);

%% Post-process

% Quantify the symmetry breaking of Mitos
% Figure 7
Nmt_front = zeros(Nt2,1);
Nmt_back  = zeros(Nt2,1);
Dir_migration = sign(Xs(Nt2,1));
for it = 1:Nt2
    xmt = Xmt(it,:);
    ymt = Ymt(it,:);
    band = abs(ymt-Xs(it,2)) <= width;
    front = (xmt - Xs(it,1))*Dir_migration >= 0;
    back = (xmt - Xs(it,1))*Dir_migration <= 0;
    Nmt_front(it) = nnz(front & band);
    Nmt_back(it)  = nnz(back & band);
end


MSDs = ((Xs(:,1)-Xs(1,1)).^2+(Xs(:,2)-Xs(1,2)).^2);

MSDmt = zeros(Nt2,Nmt);
for imt = 1:Nmt
    MSDmt(:,imt) = ((Xmt(:,imt)-Xmt(1,imt)).^2 + (Ymt(:,imt)-Ymt(1,imt)).^2); 
end
MSDMT = mean(MSDmt,2);


%% Plot

figure(3) % spindle trajectory
x_traj = Xs(:,1)*1d6;
y_traj = Xs(:,2)*1d6;
plot(x_traj,y_traj,'-*k'); hold off
set(gca, 'fontsize',13);
% axis([-20 20 -15 15]);
xlabel('x ({\mu}m)', 'fontsize',13);
ylabel('y ({\mu}m)', 'fontsize',13);
% axis equal

figure(4) % spindle and averged mitos mean square displacement
plot(Time(1:SaveTime:Nt)/60,MSDs*1d12,'-k','linewidth',2); hold on
plot(Time(1:SaveTime:Nt)/60,MSDMT*1d12,'-b','linewidth',2); hold off
set(gca, 'fontsize',13);
axis([0 round(dt*Nt/60) 0 700.1]);
xlabel('t (min)', 'fontsize',13);
ylabel('r^2 ({\mu}m^2)', 'fontsize',13);
legend('Spindle','Mitos');

% Quantify the symmetry breaking of Mitos
MSDs = (Xs(:,1)-Xs(1,1)).^2;
mask = MSDs*1d12 < 8;
nnz(mask)

Mito_ratio = Nmt_back./Nmt_front;
mask = Mito_ratio < 1.3;
nnz(mask)


figure(7) % coplot of mito ratio and spindle squared displacement
yyaxis left
plot(Time(1:SaveTime:Nt)/60,MSDs*1d12,'-','linewidth',2); 
set(gca, 'fontsize',13);
% axis([0 round(dt*Nt/60) -30 600]);
xlabel('{\it t} (min)', 'fontsize',13);
ylabel('Spindle {\it x} ^2 ({\mu}m^2)', 'fontsize',13);
yyaxis right
semilogy(Time(1:SaveTime:Nt)/60,Nmt_back./Nmt_front,'-','linewidth',0.5); 
set(gca, 'fontsize',13);
% axis([0 round(dt*Nt/60) 1 100]);
xlabel('{\it t} (min)', 'fontsize',13);
ylabel('Mito. Back to Front Ratio', 'fontsize',13);


