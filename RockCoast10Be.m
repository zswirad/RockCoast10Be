% Backward geometry-based model to explore cosmogenic 10Be concentrations 
%   across an active shore platform as a function of cliff retreat and shore 
%   platform down-wearing scenarios

% Swirad, Z. M. et al. 2020. Cosmogenic exposure dating reveals limited long-term
%   variability in erosion of a rocky coastline. Nature Communcations  11: 3804. 
%   https://doi.org/10.1038/s41467-020-17611-9

% RockCoast10Be.m
% Version 2.0
% last update: 2025-07-01 (MATLAB R2023b)
% Zuzanna Swirad (zswirad@igf.edu.pl)
% Institute of Geophysics, Polish Academy of Sciences

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

% INPUT FILES
GeomagnScalarRaw = load('geomagnetics.txt'); % Geomagnetic scalar, Sgm; c1: time (yr BP), c2: Sgm: Lifton et al., 2014
ProfileRaw = load('profile.txt'); % c1: distance from the cliff (m); c2: elevation (m a.s.l.)
RSLRaw = load('sealevel.txt'); % Relative sea level, RSL; c1: time (yr BP), c2: RSL (m a.s.l.)
Be = load('measured_inheritance.txt'); % measured 10Be concentrations; c1: distance from the cliff (m); c2: concentrations after subtracting inherited contribution (atoms/g); c3: background error and inheritance error propagating in quadrature (atoms/g) % COMMENT IF NO 10BE DATA AVAILABLE

% Interpolate geomagnetic scalar and RSL to 1 year and profile to 1 m spacing
GeomagnScalar = zeros(10001,2);
GeomagnScalar(:,1) = 0:10000;
GeomagnScalar(:,2) = interp1(GeomagnScalarRaw(:,1),GeomagnScalarRaw(:,2),0:10000);
clear GeomagnScalarRaw
Profile = zeros(max(ProfileRaw(:,1))+1,2);
Profile(:,1) = 0:max(ProfileRaw(:,1));
Profile(:,2) = interp1(ProfileRaw(:,1),ProfileRaw(:,2),0:max(ProfileRaw(:,1)));
clear ProfileRaw
RSL = zeros(max(RSLRaw(:,1))+1,2);
RSL(:,1) = 0:max(RSLRaw(:,1));
RSL(:,2) = interp1(RSLRaw(:,1),RSLRaw(:,2),0:max(RSLRaw(:,1)));
clear RSLRaw

% VARIABLES
ProdRate10Be = 4.009; % 10Be production rate (at/g/yr)
TotalTime = 7000; % total time considered (yr)
Time = 0:TotalTime; % time verctor (yr)
PlatformWidth = 300; % contemporary shore platform width (m)
DistanceCliff = 0:PlatformWidth; % distance from the cliff vector (m)
TidalRange = 4.6; % tidal range (m)
HAT = 3.2; % highest astronomical tide or highest elevation of shore platform above RSL (m a.s.l.)
LAT = -2.8; % lowest astronomical tide or lowest elevation of shore platform below RSL (m a.s.l.)
ro = 1.024; % seawater density (g/cm3)
lambda = 160; % high-energy neutrons attenuation length (g/cm2): Goesse and Phillips, 2001
lambda2 = 1.3; % attenuation length for particle flux: Dunne et al., 1999
CliffHeight = 50; % (m)
Theta = deg2rad(90); % cliff inclination angle (degrees)
dPhi = deg2rad(180); % subtended azimuth angle (degrees)
CliffWidth = round(CliffHeight/tan(Theta)); % horizontal distance between cliff base and cliff top (0 if vertical)
m_coef = 2.3; % scaling coefficient: Dunne et al., 1999
io = ProdRate10Be*(m_coef+1)/(2*pi); % incidence radiation; calculated from maximal radiation Fmax from Dunne et al., 1999
ThetaPlatform = atan(CliffHeight./(CliffWidth+DistanceCliff)); % cliff inclination across the platform
TopoShield0 = 1 - (io*dPhi/(m_coef+1)*sin(ThetaPlatform).^(m_coef+1)/ProdRate10Be);
TopoShield = zeros(size(DistanceCliff,2),2);
TopoShield(:,1) = DistanceCliff;
TopoShield(:,2) = TopoShield0;
clear TopoShield0

% CLIFF RETREAT SCENARIOS & EXPOSURE AGES
SteadyRetreatRate = [0.01 0.05:0.05:0.5]; % (m/yr)
ChangeRetreatRate = 2:2:10; % at TotalTime retreat rate was X time faster/slower than PresentRetreatRate
PresentRetreatRate = 0.1; % used to calculate the rates for acceleration/deceleration scenarios (m/yr)
ScenariosNo = length(SteadyRetreatRate) + 2*length(ChangeRetreatRate);

% Matrices follow such a scheme: top left corner refers to the cliff toe
%   and/or present time; they usually have ScenariosNo/PlatformWidth+1/TotalTime+1 
%   rows/columns:
% TotalTime+1: time BP values (Time vector)
% PlatformWidth+1: distance from the cliff values (DistanceCliff vector)
% ScenariosNo: cliff retreat scenarios

% VISUALISATION
% Re-order the rows (cliff retreat scenarios) so that they follow such a scheme:
% 1. acceleration scenarios such that the rate was X times lower at TotalTime
% 2. steady retreat rate scenarios
% 3. deceleration scenarios such that the rate was X times higher at TotalTime

cc = flipud(jet(ScenariosNo)); % colour scheme for plotting (jet scheme makes acceleration scenarios red and deceleration blue)

RetreatRate = zeros(ScenariosNo,TotalTime+1); % retreat rate through time
for m = 1:length(ChangeRetreatRate) % acceleration scenarios (from the largest change)
    end_rate = PresentRetreatRate/ChangeRetreatRate(m);
    RetreatRate(length(ChangeRetreatRate)+1-m,:) = linspace(PresentRetreatRate,end_rate,TotalTime+1);
end
for m = 1:length(SteadyRetreatRate) % steady retreat scenarios (from the fastest rate)
    RetreatRate(length(ChangeRetreatRate)+length(SteadyRetreatRate)+1-m,:) = SteadyRetreatRate(m);
end
for m = 1:length(ChangeRetreatRate) % deceleration scenarios (from the smallest change)
    end_rate = PresentRetreatRate*ChangeRetreatRate(m);
    RetreatRate(m+length(ChangeRetreatRate)+length(SteadyRetreatRate),:) = linspace(PresentRetreatRate,end_rate,TotalTime+1);
end

CliffPosition = zeros(ScenariosNo,TotalTime+1); % relative cliff position through time
for m = 1:ScenariosNo
    for n = 1:TotalTime
        CliffPosition(m,n+1) = CliffPosition(m,n) + RetreatRate(m,n);
    end
end
CliffPosition = round(CliffPosition,4);

ExpoAges = zeros(ScenariosNo,PlatformWidth+1); % cross-shore exposure ages
for m = 1:ScenariosNo
    CliffPositionScenario = CliffPosition(m,:);
    for w = 1:PlatformWidth+1
        if CliffPositionScenario(end) < (w-1)
            ExpoAges(m,w) = 0;
        else
            ExpoAges(m,w) = find(CliffPositionScenario >= (w-1),1);
        end
    end
end

figure(1)
for n = 1:ScenariosNo
    plot(Time,RetreatRate(n,:),'color',cc(n,:))
    hold on
end
hold off
title('Past cliff retreat rates')
xlabel('Time (yr BP)') 
ylabel('Cliff retreat rate (m/yr)') 

figure(2)
for m = 1:ScenariosNo
    for w = 1:PlatformWidth+1
        if ExpoAges(m,w) == 0
            ExpoAges(m,w) = ExpoAges(m,w)/0;
        end
    end
end
for n=1:ScenariosNo
    plot(DistanceCliff,ExpoAges(n,:),'color',cc(n,:))
    hold on
end
hold off
axis([0 PlatformWidth 0 TotalTime])
title('Cross-shore distribution of exposure ages')
xlabel('Distance from the cliff (m)') 
ylabel('Time (yr BP)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GEOMAGNETIC SCALAR (Sgm)
Sgm = zeros(ScenariosNo,PlatformWidth+1); % cumulative Sgm
for m = 1:ScenariosNo
    for w = 1:PlatformWidth+1
        if ExpoAges(m,w) >-1
            Sgm(m,w) = sum(GeomagnScalar(1:ExpoAges(m,w),2))/ExpoAges(m,w);
        else
            Sgm(m,w) = nan;
        end
    end
end

figure(3)
plot(GeomagnScalar(:,1),GeomagnScalar(:,2),'k')
title('Geomagnetic scalar (Lifton et al., 2014)')
xlabel('Time (yr BP)') 
ylabel('Geomagnetic scalar, Sgm')

figure(4)
for n=1:ScenariosNo
    plot(DistanceCliff,Sgm(n,:),'color',cc(n,:))
    hold on
end
hold off
axis([0 PlatformWidth 0.9 1.15]) 
title('Cross-shore distribution of total geomagnetic scalar')
xlabel('Distance from the cliff (m)') 
ylabel('Geomagnetic scalar, Sgm') 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TOPOGRAPHIC SHIELDING (Stopo)
Stopo = zeros(ScenariosNo,PlatformWidth+1); % cumulative Stopo
for m = 1:ScenariosNo
    DistanceCliffScenario = zeros(TotalTime+1,PlatformWidth+1);
    StopoScenario = zeros(TotalTime+1,PlatformWidth+1);
    for w = 1:PlatformWidth+1
        for n = 1:TotalTime+1
            if ExpoAges(m,w) >= n
                DistanceCliffScenario(n,w) = DistanceCliff(w)-CliffPosition(m,n);
                if DistanceCliffScenario(n,w) <= 0
                    StopoScenario(n,w) = TopoShield(1,2);
                else
                    for w2 = 2:PlatformWidth+1
                        if DistanceCliffScenario(n,w) > TopoShield(w2-1,1) && DistanceCliffScenario(n,w) <= TopoShield(w2,1)
                            StopoScenario(n,w) = (DistanceCliffScenario(n,w)-TopoShield(w2-1,1)) * (TopoShield(w2,2)-TopoShield(w2-1,2)) + TopoShield(w2-1,2);
                        end
                    end
                end
            end
        end
    end
    Stopo(m,:) = sum(StopoScenario,1)./ExpoAges(m,:);
end

figure(5)
plot(TopoShield(:,1),TopoShield(:,2),'k')
title('Present topographic shielding')
xlabel('Distance from the cliff (m)') 
ylabel('Topographic shielding, Stopo')

figure(6)
for n = 1:ScenariosNo
    plot(DistanceCliff,Stopo(n,:),'color',cc(n,:))
    hold on
end
hold off
axis([0 PlatformWidth 0.5 1])
title('Cross-shore distribution of cumulative topographic shielding')
xlabel('Distance from the cliff (m)') 
ylabel('Topographic shielding, Stopo') 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SHORE PLATFORM EROSION SCALAR (Ser) && WATER SHIELDING (Sw)
% Assumption: shore platform cannot be > HAT or < LAT relative to the sea level 
% Platform erosion scenarios:
% 1. zero platform erosion
% 2. steady-state model
% 3. platform widening model

BestFit = polyfit(Profile(:,1),Profile(:,2),1);
BestFitProfile = BestFit(1)*DistanceCliff + BestFit(2);

figure(7)
plot(Profile(:,1),Profile(:,2),'k')
hold on 
plot(Profile(:,1),BestFitProfile,'r')
hold on 
yline(HAT,'--k')
hold on 
yline(LAT,'--r')
hold off
legend('real profile','best-fit profile','HAT','LAT')
title('Cross-shore profile')
xlabel('Distance from the cliff (m)') 
ylabel('Elevation (m a.s.l.)') 

figure(8)
plot(RSL(:,1),RSL(:,2),'k')
title('Relative sea level')
xlabel('Time (yr BP)') 
ylabel('Elevation (m a.s.l.)')

tidal_time = 0:0.01:24;
tidal_level = 0.5 * TidalRange * sin(tidal_time*2.*pi/12.);
Tides = zeros(round((HAT-LAT)*10)+1,3); % c1: elevation (m a.s.l.), c2: tidal duration distribution (adds to 1), c3: water shielding
TidesBin = 0.1; % (m)
Tides(:,1) = HAT:-TidesBin:LAT;
for n = 1:size(tidal_level,2)
    for m = 1:size(Tides,1)
        if tidal_level(n) > Tides(m,1)-TidesBin && tidal_level(n) <= Tides(m,1)
            Tides(m,2) = Tides(m,2) + 1;
        end
    end
end
Tides(:,2) = Tides(:,2)./size(tidal_level,2);

for m = 1:size(Tides,1)
    Sw_temp = zeros(size(Tides,1),1);
    for n = 1:size(Tides,1)
        if Tides(m,1) >= Tides(n,1)
            Sw_temp(n) = 1;
        else
            WaterDepth = (Tides(n,1) - Tides(m,1)) * 100; % (cm)
            Sw_temp(n) = exp(-ro*WaterDepth/lambda);
        end
    end
    Sw_temp = Sw_temp.*Tides(:,2);
    Tides(m,3) = sum(Sw_temp);
end

figure(9)
plot(Tides(:,3),Tides(:,1),'k')
title('Present water shielding')
xlabel('Water shielding, Sw') 
ylabel('Elevation (m a.s.l.)')
axis([0 1 LAT HAT])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Scenario 1. zero platform erosion (output: Sw1; Ser1 = 1)

Sw1 = zeros(ScenariosNo,PlatformWidth+1); % cumulative Sw1
for m = 1:ScenariosNo
    Sw1Scenario = zeros(TotalTime+1,PlatformWidth+1);
    for n = 1:TotalTime+1
        ElevRSL = BestFitProfile' - RSL(n,2); % elevation relative to sea level at time n
        for w = 1:PlatformWidth+1
            if ElevRSL(w) > RSL(n,2) + HAT || ElevRSL(w) < RSL(n,2) + LAT
                ElevRSL(w) = nan;
            end
            if w-1 >= CliffPosition(m,n) 
                if isnan(ElevRSL(w)) == 0
                    for t = 1:size(Tides,1)
                        if ElevRSL(w) > Tides(t,1)-TidesBin && ElevRSL(w) <= Tides(t,1)
                            Sw1Scenario(n,w) = Tides(t,3);
                        end
                    end
                else
                    Sw1Scenario(n,w) = nan;
                end
            end
        end
    end
    for w = 1:PlatformWidth+1
        if any(isnan(Sw1Scenario(:,w))) == 0
            Sw1(m,w) = (sum(Sw1Scenario(:,w))) /ExpoAges(m,w);
        else
            Sw1(m,w) = nan;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Scenario 2. steady-state model (outputs: Sw2, Ser2)

Elev2_noRSL = zeros(ScenariosNo,TotalTime+1); % elevation relative to the present without RSL
Down2_noRSL = RetreatRate*atan((BestFitProfile(1)-BestFitProfile(end))/PlatformWidth); % shore platform down-wearing rate without RSL
for m = 1:ScenariosNo
    for n = 2:TotalTime+1
        Elev2_noRSL(m,n) = Elev2_noRSL(m,n-1) + Down2_noRSL(m,n-1);
    end
end

Sw2 = zeros(ScenariosNo,PlatformWidth+1); % cumulative Sw2
Elev2_RSL = zeros(ScenariosNo,TotalTime+1); % elevation relative to the present including RSL
for m = 1:ScenariosNo
    Sw2Scenario = zeros(TotalTime+1,PlatformWidth+1);
    Elev2_RSL(m,:) = Elev2_noRSL(m,:) + RSL(:,2)';
    for n = 1:TotalTime+1
        Elev2_RSL_profile = BestFitProfile + Elev2_RSL(m,n);
        for w = 1:PlatformWidth+1
            if Elev2_RSL_profile(w) > RSL(n,2) + HAT || Elev2_RSL_profile(w) < RSL(n,2) + LAT
                Elev2_RSL_profile(w) = nan;
            end
            if w-1 >= CliffPosition(m,n)
                if isnan(Elev2_RSL_profile(w)) == 0
                    for t = 1:size(Tides,1)
                        if Elev2_RSL_profile(w) > Tides(t,1)-TidesBin && Elev2_RSL_profile(w) <= Tides(t,1)
                            Sw2Scenario(n,w) = Tides(t,3);
                        end
                    end
                else
                    Sw2Scenario(n,w) = nan;
                end
            end
        end
    end
    for w = 1:PlatformWidth+1
        if any(isnan(Sw2Scenario(:,w))) == 0
            Sw2(m,w) = (sum(Sw2Scenario(:,w))) /ExpoAges(m,w);
        else
            Sw2(m,w) = nan;
        end
    end
 end

Ser2Scenario = exp(-Elev2_RSL/lambda2);
Ser2Scenario(Ser2Scenario>1) = 1;

Ser2 = zeros(ScenariosNo,PlatformWidth+1); % cumulative Ser2
for m = 1:ScenariosNo
    for n = 1:TotalTime+1
        for w = 1:PlatformWidth+1
            if isnan(Sw2(m,w)) == 0
                Ser2(m,w) = sum(Ser2Scenario(m,1:ExpoAges(m,w)))/ExpoAges(m,w);
            else
                Ser2(m,w) = nan;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Scenario 3. platform widening model (outputs: Sw3, Ser3)
SlopeTemporal = zeros(ScenariosNo,TotalTime+1); % slope through time
PlatformWidthTemporal = PlatformWidth - CliffPosition; % platform width through time
for m = 1:ScenariosNo
    for n = 1:TotalTime+1
        if PlatformWidthTemporal(m,n) < 0
            PlatformWidthTemporal(m,n) = 0;
        else
            SlopeTemporal(m,n) = atan((BestFitProfile(1)-BestFitProfile(end))/PlatformWidthTemporal(m,n)); % in radians
        end
    end
end

Ser3 = zeros(ScenariosNo,PlatformWidth+1); % cumulative Ser3
Sw3 = zeros(ScenariosNo,PlatformWidth+1); % cumulative Sw3

for m = 1:ScenariosNo
    Elev3_Abs = zeros(TotalTime+1,PlatformWidth+1); % elevation relative to the present through time
    Elev3_Abs(1,:) = BestFitProfile;
    Elev3_RSL = zeros(TotalTime+1,PlatformWidth+1); % elevation relative to the RSL through time
    Elev3_RSL(1,:) = BestFitProfile;
    Sw3Scenario = zeros(TotalTime+1,PlatformWidth+1);
    SampleDepth3 = zeros(TotalTime+1,PlatformWidth+1); % depth of presently exposed surface through time
    Ser3Scenario = zeros(TotalTime+1,PlatformWidth+1);
    Ser3Scenario(1,:) = 1;
    for n = 2:TotalTime+1
        for w = 1:PlatformWidth+1
            if w-1 >= CliffPosition(m,n) && SlopeTemporal(m,n) ~= 0
                Elev3_Abs(n,w) = tan(SlopeTemporal(m,n)) * (PlatformWidth+1-w) + BestFitProfile(end); 
                SampleDepth3(n,w) = Elev3_Abs(n,w) - BestFitProfile(w);
                Ser3Scenario(n,w) = exp(-SampleDepth3(n,w) /lambda2);
                Elev3_RSL(n,w) = Elev3_Abs(n,w) + RSL(n,2);
                if Elev3_RSL(n,w) > RSL(n,2) + HAT || Elev3_RSL(n,w) < RSL(n,2) + LAT
                    Elev3_RSL(n,w) = nan;
                end
                if w-1 >= CliffPosition(m,n)
                    if isnan(Elev3_RSL(n,w)) == 0
                        for t = 1:size(Tides,1)
                            if Elev3_RSL(n,w) > Tides(t,1)-TidesBin && Elev3_RSL(n,w) <= Tides(t,1)
                                Sw3Scenario(n,w) = Tides(t,3);
                            end
                        end
                    else
                        Sw3Scenario(n,w) = nan;
                    end
                end
            end
        end
    end
    for w = 1:PlatformWidth+1
        if any(isnan(Sw3Scenario(:,w))) == 0
            Ser3(m,w) = (sum(Ser3Scenario(:,w)))/ExpoAges(m,w);
            Sw3(m,w) = (sum(Sw3Scenario(:,w)))/ExpoAges(m,w);
        else
            Ser3(m,w) = nan;
            Sw2(m,w) = nan;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(10)

subplot(1,3,1)
for n=1:ScenariosNo
    plot(DistanceCliff,Sw1(n,:),'color',cc(n,:))
    hold on
end
hold off
axis([0 PlatformWidth 0 1])
title('Zero platform erosion')
xlabel('Distance from the cliff (m)') 
ylabel('Water shielding, Sw') 

subplot(1,3,2)
for n = 1:ScenariosNo
    plot(DistanceCliff,Sw2(n,:),'color',cc(n,:))
    hold on
end
hold off
axis([0 PlatformWidth 0 1])
title('Steady-state model')
xlabel('Distance from the cliff (m)') 
ylabel('Water shielding, Sw') 

subplot(1,3,3)
for n=1:ScenariosNo
    plot(DistanceCliff,Sw3(n,:),'color',cc(n,:))
    hold on
end
hold off
axis([0 PlatformWidth 0 1])
title('Platform widening model')
xlabel('Distance from the cliff (m)') 
ylabel('Water shielding, Sw') 

figure(11)

subplot(1,2,1)
for n = 1:ScenariosNo
    plot(DistanceCliff,Ser2(n,:),'color',cc(n,:))
    hold on
end
hold off
axis([0 PlatformWidth 0 1])
title('Steady-state model')
xlabel('Distance from the cliff (m)') 
ylabel('Platform erosion scalar, Ser') 

subplot(1,2,2)
for n=1:ScenariosNo
    plot(DistanceCliff,Ser3(n,:),'color',cc(n,:))
    hold on
end
hold off
axis([0 PlatformWidth 0 1]) 
title('Platform widening model')
xlabel('Distance from the cliff (m)') 
ylabel('Platform erosion scalar, Ser') 

% Calculate modelled concentrations
conc1 = ExpoAges.*Stopo.*Sgm.*Sw1.*ProdRate10Be; % zero platform erosion
conc2 = ExpoAges.*Stopo.*Sgm.*Sw2.*Ser2.*ProdRate10Be; % steady-staten model
conc3 = ExpoAges.*Stopo.*Sgm.*Sw3.*Ser3*ProdRate10Be; % platform widening model

%%% COMMENT TILL THE END IF NO 10BE DATA AVAILABLE

% Visualise raw concentrations and topography
figure(12)
yyaxis left
errorbar(Be(:,1),Be(:,2),Be(:,3),'color','k','LineStyle','none');
hold on
scatter(Be(:,1),Be(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0])
yyaxis right
plot(DistanceCliff,Profile(:,2),'-r')
yyaxis left
title('Cross-shore distribution of measured inheritance-corrected 10Be concentrations')
xlabel('Distance from the cliff (m)') 
ylabel('10Be concentrations (atoms/g)') 
yyaxis right
ylabel('Elevation (m a.s.l.)')

% Visualise measured 10Be concentrations and those modelled with all scenarios  
figure(13)
subplot(1,3,1)
for n=1:ScenariosNo
    plot(DistanceCliff,conc1(n,:),'color',cc(n,:))
    hold on
end
errorbar(Be(:,1),Be(:,2),Be(:,3),'color','k','LineStyle','none');
hold on
scatter(Be(:,1),Be(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0])
hold off
axis([0 PlatformWidth 0 max(Be(:,2))+2000])
title('Zero platform erosion')
xlabel('Distance from the cliff (m)') 
ylabel('10Be concentrations (atoms/g)') 

subplot(1,3,2)
for n=1:ScenariosNo
    plot(DistanceCliff,conc2(n,:),'color',cc(n,:))
    hold on
end
errorbar(Be(:,1),Be(:,2),Be(:,3),'color','k','LineStyle','none');
hold on
scatter(Be(:,1),Be(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0])
hold off
axis([0 PlatformWidth 0 max(Be(:,2))+2000])
title('Steady-state model')
xlabel('Distance from the cliff (m)') 
ylabel('10Be concentrations (atoms/g)') 

subplot(1,3,3)
for n=1:ScenariosNo
    plot(DistanceCliff,conc3(n,:),'color',cc(n,:))
    hold on
end
errorbar(Be(:,1),Be(:,2),Be(:,3),'color','k','LineStyle','none');
hold on
scatter(Be(:,1),Be(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0])
hold off
axis([0 PlatformWidth 0 max(Be(:,2))+2000])
title('Platform widening model')
xlabel('Distance from the cliff (m)') 
ylabel('10Be concentrations (atoms/g)') 

% Find the best-fit scenario using mean squared difference
conc1_be = conc1(:,Be(:,1));
conc2_be = conc2(:,Be(:,1));
conc3_be = conc3(:,Be(:,1));

for n = 1: size(conc1_be,1)
    conc1_msd(n,:) = (Be(:,2) - conc1_be(n,:)').^2;
    conc2_msd(n,:) = (Be(:,2) - conc2_be(n,:)').^2;
    conc3_msd(n,:) = (Be(:,2) - conc3_be(n,:)').^2;
end

figure(14)
subplot(1,3,1)
for n=1:ScenariosNo
    plot(Be(:,1),conc1_msd(n,:),'color',cc(n,:))
    hold on
end
hold off
xlim([0 PlatformWidth])
title('Zero platform erosion')
xlabel('Distance from the cliff (m)') 
ylabel('Squared difference') 

subplot(1,3,2)
for n=1:ScenariosNo
    plot(Be(:,1),conc2_msd(n,:),'color',cc(n,:))
    hold on
end
hold off
xlim([0 PlatformWidth])
title('Steady-state model')
xlabel('Distance from the cliff (m)') 
ylabel('Squared difference') 

subplot(1,3,3)
for n=1:ScenariosNo
    plot(Be(:,1),conc3_msd(n,:),'color',cc(n,:))
    hold on
end
hold off
xlim([0 PlatformWidth])
title('Platform widening model')
xlabel('Distance from the cliff (m)') 
ylabel('Squared difference') 

mse1 = mean(conc1_msd,2);
mse2 = mean(conc2_msd,2);
mse3 = mean(conc3_msd,2);

figure(15)
for n=1:ScenariosNo
    scatter(n,mse1(n),[],cc(n,:),'filled')
    scatter(n,mse2(n),[],cc(n,:))
    scatter(n,mse3(n),[],cc(n,:),'*')
    hold on
end
xline(size(ChangeRetreatRate,2)+0.5)
xline(ScenariosNo-size(ChangeRetreatRate,2)+0.5)
hold off
xlim([1 ScenariosNo])
title({'Mean squared difference','(filled: zero platform erosion, o: steady-state, *: platform widening)'})
xlabel('Scenario (acceleration -> steady -> deceleration)') 
ylabel('Mean squared difference') 

mse_full = [mse1 mse2 mse3];
min_value = min(min(mse_full));
[x,y]=find(mse_full==min_value);

disp('The best-fit scenario:')
if y==1
    disp('zero platform erosion')
elseif y==2
    disp('steady-state model')
else
    disp('platform widening model')
end
disp('present cliff retreat rate (m/yr):')
RetreatRate(x,1)
disp('initial cliff retreat rate (m/yr):')
RetreatRate(x,end)
if RetreatRate(x,1) == RetreatRate(x,end)
    disp('steady cliff retreat')
elseif RetreatRate(x,1) > RetreatRate(x,end)
    disp('acceleration')
else
    disp('deceleration')
end

figure(16)
if y==1
    plot(DistanceCliff,conc1(x,:),'k')
elseif y==2
    plot(DistanceCliff,conc2(x,:),'k')
else
    plot(DistanceCliff,conc3(x,:),'k')
end
hold on
errorbar(Be(:,1),Be(:,2),Be(:,3),'color','k','LineStyle','none');
hold on
scatter(Be(:,1),Be(:,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0])
hold off
xlim([0 PlatformWidth])
title('Cross-shore distribution of measured and modelled 10Be concentrations')
xlabel('Distance from the cliff (m)') 
ylabel('10Be concentrations (atoms/g)') 