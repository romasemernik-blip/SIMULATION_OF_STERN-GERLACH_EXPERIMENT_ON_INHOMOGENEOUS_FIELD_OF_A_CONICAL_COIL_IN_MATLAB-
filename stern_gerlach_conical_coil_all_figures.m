function stern_gerlach_conical_coil_all_figures()
% ЕДИНАЯ МОДЕЛЬ ЭКСПЕРИМЕНТА ШТЕРНА-ГЕРЛАХА С КОНИЧЕСКОЙ КАТУШКОЙ
% Версия согласована с описанием разделов 3.1–5 реферата.
%
% Геометрия: Y — ось катушки и направление атомного пучка;
% X — горизонтальная координата; Z — вертикальная координата детектора.
%
% ПОЛЕ ПО МЕТОДИКЕ РЕФЕРАТА (раздел 3.2):
% 1. Поле обмотки — сумма вкладов кольцевых слоёв.
% 2. Радиус кольцевого слоя линейно меняется вдоль Y.
% 3. Осевая компонента отдельного кольцевого витка:
%    dBy = mu0*n*I*ds*R^2/(2*(R^2+(Y-s)^2)^(3/2)).
% 4. dBy/dY рассчитывается конечными разностями.
% 5. Поперечная компонента в параксиальном приближении:
%    Br = -(r/2)*dBy/dY.
% 6. |B| = sqrt(By^2 + Br^2), затем сила для двух проекций:
%    F_s = s*muB*grad(|B|).
%
% Это именно приближённая осесимметричная схема из реферата, а НЕ точное
% поле каждой винтовой нити обмотки. Не смешивать её с моделью эллиптических
% интегралов для точного поля круговых токовых колец.
%
% ВАЖНЫЕ ОГРАНИЧЕНИЯ:
% - атомы независимы; столкновения и рассеяние на остаточном газе не заданы;
% - искусственные смещения и импульсы после катушки отсутствуют;
% - модель двух проекций магнитного момента — адиабатическое приближение;
% - если расчёт не даёт заметного разделения, программа его не дорисовывает.
%
% Запуск: stern_gerlach_conical_coil_all_figures
% Все рисунки, CSV и MAT сохраняются в папку stern_gerlach_figures.

close all; clc;
rng(7,'twister');

%% 1. ПАРАМЕТРЫ КОНИЧЕСКОЙ КАТУШКИ (таблица 3.1.1 реферата)
p.mu0 = 4*pi*1e-7;       % Гн/м
p.n = 3200;              % плотность витков, 1/м
p.I = 650;              % ток, А
p.a = 0.32;              % радиус входного сечения, м
p.b = 0.44;              % радиус выходного сечения, м
p.y0 = 0.00;             % координата начала катушки, м
p.L = 1.15;              % длина катушки, м
p.y1 = p.y0+p.L;
p.Nlayers = 720;         % узлы квадратуры по длине, не физическое число витков
p.NphysicalTurns = round(p.n*p.L);

%% 2. ОБЛАСТЬ РАСЧЁТА ПОЛЯ
rGrid = linspace(0,0.55,241);
yGrid = linspace(-0.35,2.30,701);
[RR,YY] = ndgrid(rGrid,yGrid);
turnY = linspace(p.y0,p.y1,p.Nlayers);
turnR = p.a+(p.b-p.a)*(turnY-p.y0)/p.L;
ds = (p.y1-p.y0)/(p.Nlayers-1);
weights = ds*ones(size(turnY));
weights([1 end]) = ds/2;

fprintf('Расчёт поля по методике раздела 3.2 реферата...\n');
By = zeros(size(RR));
% Суммируем осевое поле кольцевых слоёв по формуле кругового витка.
for k = 1:p.Nlayers
    dyLocal = YY-turnY(k);
    By = By + p.mu0*p.n*p.I*weights(k)*turnR(k)^2 ./ ...
        (2*(turnR(k)^2+dyLocal.^2).^(3/2));
end

% Центральные конечные разности по Y; на границах — односторонние.
dBydY = zeros(size(By));
dY = yGrid(2)-yGrid(1);
dBydY(:,2:end-1) = (By(:,3:end)-By(:,1:end-2))/(2*dY);
dBydY(:,1) = (By(:,2)-By(:,1))/dY;
dBydY(:,end) = (By(:,end)-By(:,end-1))/dY;

% Параксиальная поперечная компонента из формулы раздела 3.2.
Br = -0.5.*RR.*dBydY;
Bmag = hypot(By,Br);

% Градиент модуля поля для расчёта силы на атомы.
[dBdr,dBdy] = gradient(Bmag,rGrid(2)-rGrid(1),dY);

%% 3. АТОМЫ СЕРЕБРА И ИСХОДНЫЙ ПУЧОК (раздел 4.1)
N_atoms = 100000;
mAg = 1.791e-25;             % кг
muB = 9.2740100783e-24;      % Дж/Тл
kB = 1.380649e-23;           % Дж/К
T = 1000;                    % К
sourceY = -0.34;
sigmaSource = 0.25e-3;       % малый нормальный разброс координат
coneHalfAngleDeg = 1.0;      % направленный пучок вдоль +Y

% Характерная скорость Максвелла—Больцмана.
vMeanMB = sqrt(8*kB*T/(pi*mAg));
% В реферате дополнительно указан приближённый диапазон 435–525 м/с.
% Используем этот диапазон как модельный разброс вокруг характерной скорости.
vMin = 435; vMax = 525;
speed0 = vMin+(vMax-vMin)*rand(N_atoms,1);

% Независимое направление внутри малого телесного угла вокруг +Y.
cosTheta = 1-rand(N_atoms,1)*(1-cosd(coneHalfAngleDeg));
sinTheta = sqrt(max(0,1-cosTheta.^2));
phi = 2*pi*rand(N_atoms,1);
vx = speed0.*sinTheta.*cos(phi);
vz = speed0.*sinTheta.*sin(phi);
vy = speed0.*cosTheta;

x = sigmaSource*randn(N_atoms,1);
z = sigmaSource*randn(N_atoms,1);
y = sourceY*ones(N_atoms,1);

% Две равновероятные проекции магнитного момента.
spin = ones(N_atoms,1);
spin(rand(N_atoms,1)<0.5) = -1;

%% 4. ДВЕ ЩЕЛИ И ТРИ ДЕТЕКТОРА (раздел 4.1)
slit1Y = -0.30;
slit1Radius = 1.5e-3;        % первая щель — окружность (формула 13)
slit2Y = -0.10;
slit2HalfWidth = 0.75e-3;
slit2HalfHeight = 2.5e-3;
slit2Zshift = 0;             % пересечение двух парабол (формула 14)
detectorY = [1.46 1.80 2.20];
detectorHalfX = 0.08;
detectorHalfZ = 0.08;

%% 5. ИНТЕРПОЛЯТОРЫ ПОЛЯ И ГРАДИЕНТА
% Матрицы имеют размер length(rGrid) x length(yGrid).
% interp2 ожидает сетку Y по горизонтали и r по вертикали.
F_r = @(rq,yq) interp2(yGrid,rGrid,dBdr,yq,rq,'linear');
F_y = @(rq,yq) interp2(yGrid,rGrid,dBdy,yq,rq,'linear');

%% 6. РИСУНКИ МАГНИТНОГО ПОЛЯ — РАЗДЕЛ 3.2
outDir = 'stern_gerlach_figures';
if ~exist(outDir,'dir'), mkdir(outDir); end

% Рисунок 1: трёхмерная геометрия конической катушки.
fig = figure('Color','w','Name','3D геометрия конической катушки');
hold on; axis equal; grid on; view(35,22);
ang = linspace(0,2*pi,120);
for q = 1:45
    yy = p.y0+(p.y1-p.y0)*(q-1)/44;
    rq = p.a+(p.b-p.a)*(yy-p.y0)/p.L;
    plot3(rq*cos(ang),yy*ones(size(ang)),rq*sin(ang),'Color',[0.1 0.35 0.75]);
end
plot3([0 0],[p.y0 p.y1],[0 0],'k--','LineWidth',1.2);
xlabel('X, м'); ylabel('Y, м — ось пучка'); zlabel('Z, м');
title('Коническая катушка, соосная атомному пучку');
exportgraphics(fig,fullfile(outDir,'01_coil_3d.png'),'Resolution',220);

% Рисунок 2: продольное сечение конической обмотки.
fig = figure('Color','w','Name','Продольное сечение катушки');
hold on; axis equal; grid on;
yy = linspace(p.y0,p.y1,60);
rr = p.a+(p.b-p.a)*(yy-p.y0)/p.L;
for q=1:numel(yy)
    plot([-rr(q) rr(q)],[yy(q) yy(q)],'Color',[0.15 0.4 0.8]);
end
plot(rr, yy, 'k--');
plot(-rr, yy, 'k--');
xlabel('Радиальная координата r, м'); ylabel('Y, м');
title('Сечение конической обмотки: радиус витков меняется вдоль Y');
exportgraphics(fig,fullfile(outDir,'02_coil_cross_section.png'),'Resolution',220);

% Рисунок 3: карта модуля поля и осевой компоненты.
fig = figure('Color','w','Name','Распределение магнитного поля');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile;
contourf(yGrid,rGrid,Bmag,36,'LineColor','none'); colorbar;
hold on; plot([p.y0 p.y1],[p.a p.b],'w--','LineWidth',1.3);
xlabel('Y, м'); ylabel('r, м'); title('|B(r,Y)|, Тл');
nexttile;
contourf(yGrid,rGrid,By,36,'LineColor','none'); colorbar;
hold on; plot([p.y0 p.y1],[p.a p.b],'w--','LineWidth',1.3);
xlabel('Y, м'); ylabel('r, м'); title('B_y(r,Y), Тл');
exportgraphics(fig,fullfile(outDir,'03_field_maps.png'),'Resolution',220);

% Рисунок 4: поперечная компонента и осевой градиент.
fig = figure('Color','w','Name','Компоненты и градиент поля');
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile; contourf(yGrid,rGrid,Br,36,'LineColor','none'); colorbar;
xlabel('Y, м'); ylabel('r, м'); title('B_r(r,Y), Тл');
nexttile; contourf(yGrid,rGrid,dBdy,36,'LineColor','none'); colorbar;
xlabel('Y, м'); ylabel('r, м'); title('∂|B|/∂Y, Тл/м');
Baxis = By(1,:);
dBaxis = dBydY(1,:);
nexttile; plot(yGrid,Baxis,'LineWidth',1.6); grid on;
xline(p.y0,'r--'); xline(p.y1,'r--');
xlabel('Y, м'); ylabel('B_y(0,Y), Тл'); title('Осевое поле');
nexttile; plot(yGrid,dBaxis,'LineWidth',1.6); grid on;
xline(p.y0,'r--'); xline(p.y1,'r--');
xlabel('Y, м'); ylabel('dB_y/dY, Тл/м'); title('Градиент осевой компоненты');
exportgraphics(fig,fullfile(outDir,'04_field_components_and_gradient.png'),'Resolution',220);

%% 7. ОБЩАЯ СХЕМА УСТАНОВКИ
fig = figure('Color','w','Name','Общая схема установки');
hold on; grid on; axis equal; view(35,20);
[th,ss] = meshgrid(linspace(0,2*pi,60),linspace(0,1,30));
rad = p.a+(p.b-p.a)*ss;
surf(rad.*cos(th),p.y0+p.L*ss,rad.*sin(th),'FaceAlpha',0.10,'EdgeAlpha',0.10);
plot3([0 0],[sourceY detectorY(end)],[0 0],'k--');
scatter3(0,sourceY,0,45,'filled'); text(0,sourceY,0.04,'Источник Ag');
for yy0=[slit1Y slit2Y]
    plot3([-0.012 0.012],[yy0 yy0],[0 0],'k','LineWidth',3);
end
text(0,slit1Y,0.04,'Щель 1'); text(0,slit2Y,0.04,'Щель 2');
for q=1:3
    plot3([-detectorHalfX detectorHalfX],[detectorY(q) detectorY(q)],[0 0],'r','LineWidth',2);
    text(0,detectorY(q),0.04,sprintf('Детектор %d',q));
end
xlabel('X, м'); ylabel('Y, м'); zlabel('Z, м');
title('Источник, две щели, катушка и три детектора');
exportgraphics(fig,fullfile(outDir,'05_experimental_setup.png'),'Resolution',220);

%% 8. ДИНАМИКА АТОМОВ — ЯВНЫЙ МЕТОД ЭЙЛЕРА (раздел 4.2)
fprintf('Моделирование %d атомов серебра...\n',N_atoms);
fprintf('Средняя скорость Максвелла при T=%g K: %.2f м/с; модельный диапазон: %.0f–%.0f м/с\n', ...
    T,vMeanMB,vMin,vMax);
fprintf('Физическое число витков nL ≈ %d; узлов квадратуры: %d\n', ...
    p.NphysicalTurns,p.Nlayers);

alive = true(N_atoms,1);
slit1Done = false(N_atoms,1); slit2Done = false(N_atoms,1);
passed1 = false(N_atoms,1); passed2 = false(N_atoms,1);
hitX = cell(3,1); hitZ = cell(3,1);
for q=1:3
    hitX{q}=nan(N_atoms,1); hitZ{q}=nan(N_atoms,1);
end

dt = 1e-5;  % с — как указано в реферате
maxSteps = ceil((detectorY(end)-sourceY)/(max(100,vMeanMB*0.5)*dt)*1.4);
maxSteps = min(maxSteps,15000);
trajCount = 40;
sampleIds = randperm(N_atoms,min(trajCount,N_atoms));
trajX=nan(maxSteps+1,numel(sampleIds));
trajY=trajX; trajZ=trajX;
trajX(1,:)=x(sampleIds); trajY(1,:)=y(sampleIds); trajZ(1,:)=z(sampleIds);

for it=1:maxSteps
    ids=find(alive);
    if isempty(ids), break; end
    xOld=x(ids); yOld=y(ids); zOld=z(ids);
    rOld=hypot(xOld,zOld);
    rq=min(max(rOld,rGrid(1)),rGrid(end));
    yq=min(max(yOld,yGrid(1)),yGrid(end));
    gr=F_r(rq,yq); gy=F_y(rq,yq);
    gr(~isfinite(gr))=0; gy(~isfinite(gy))=0;
    ex=zeros(size(rOld)); ez=ex;
    nz=rOld>1e-14; ex(nz)=xOld(nz)./rOld(nz); ez(nz)=zOld(nz)./rOld(nz);

    % U_s=-s*mu_B*|B|; F_s=s*mu_B*grad(|B|).
    ax=spin(ids)*(muB/mAg).*gr.*ex;
    ay=spin(ids)*(muB/mAg).*gy;
    az=spin(ids)*(muB/mAg).*gr.*ez;

    % Явный Эйлер: обе координаты используют старую скорость.
    vxOld=vx(ids); vyOld=vy(ids); vzOld=vz(ids);
    xNew=xOld+vxOld*dt; yNew=yOld+vyOld*dt; zNew=zOld+vzOld*dt;
    vx(ids)=vxOld+ax*dt; vy(ids)=vyOld+ay*dt; vz(ids)=vzOld+az*dt;
    x(ids)=xNew; y(ids)=yNew; z(ids)=zNew;

    % Первая щель: окружность (13), проверка при пересечении плоскости.
    c1=~slit1Done(ids) & yOld<slit1Y & yNew>=slit1Y;
    if any(c1)
        loc=find(c1); f=(slit1Y-yOld(loc))./(yNew(loc)-yOld(loc));
        xx=xOld(loc)+f.*(xNew(loc)-xOld(loc));
        zz=zOld(loc)+f.*(zNew(loc)-zOld(loc));
        ok=xx.^2+zz.^2<=slit1Radius^2;
        g=ids(loc); passed1(g(ok))=true;
        slit1Done(g)=true; alive(g(~ok))=false;
    end

    % Вторая щель: пересечение двух парабол (14).
    c2=~slit2Done(ids) & yOld<slit2Y & yNew>=slit2Y;
    if any(c2)
        loc=find(c2); f=(slit2Y-yOld(loc))./(yNew(loc)-yOld(loc));
        xx=xOld(loc)+f.*(xNew(loc)-xOld(loc));
        zz=zOld(loc)+f.*(zNew(loc)-zOld(loc));
        hloc=slit2HalfHeight.*(1-(xx./slit2HalfWidth).^2);
        ok=abs(xx)<=slit2HalfWidth & abs(zz-slit2Zshift)<=hloc;
        g=ids(loc); passed2(g(ok))=true;
        slit2Done(g)=true; alive(g(~ok))=false;
    end

    % Три вертикальные детекторные плоскости; попадание интерполируется
    % между соседними шагами, как описано в разделе 4.2.
    for q=1:3
        pending=~isfinite(hitX{q}(ids));
        cross=pending & yOld<detectorY(q) & yNew>=detectorY(q);
        if any(cross)
            loc=find(cross); f=(detectorY(q)-yOld(loc))./(yNew(loc)-yOld(loc));
            xx=xOld(loc)+f.*(xNew(loc)-xOld(loc));
            zz=zOld(loc)+f.*(zNew(loc)-zOld(loc));
            ok=abs(xx)<=detectorHalfX & abs(zz)<=detectorHalfZ;
            g=ids(loc(ok)); hitX{q}(g)=xx(ok); hitZ{q}(g)=zz(ok);
        end
    end
    alive(ids(yNew>=detectorY(end)))=false;
    trajX(it+1,:)=x(sampleIds); trajY(it+1,:)=y(sampleIds); trajZ(it+1,:)=z(sampleIds);
end
last=find(any(isfinite(trajY),2),1,'last');
if isempty(last), last=1; end
trajX=trajX(1:last,:); trajY=trajY(1:last,:); trajZ=trajZ(1:last,:);

%% 9. РИСУНОК ТРАЕКТОРИЙ
fig=figure('Color','w','Name','Траектории атомов серебра');
hold on; grid on; view(35,22);
for q=1:size(trajX,2)
    plot3(trajX(:,q)*1e3,trajY(:,q),trajZ(:,q)*1e3,'LineWidth',0.7);
end
xlabel('X, мм'); ylabel('Y, м'); zlabel('Z, мм');
title('Траектории атомов серебра в поле конической катушки');
exportgraphics(fig,fullfile(outDir,'06_atom_trajectories_3d.png'),'Resolution',220);

%% 10. РИСУНОК РАСПРЕДЕЛЕНИЙ НА ТРЁХ ДЕТЕКТОРАХ
fig=figure('Color','w','Name','Попадания на детекторах');
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for q=1:3
    nexttile; hold on; grid on;
    valid=isfinite(hitX{q}) & isfinite(hitZ{q});
    ip=valid & spin==1; im=valid & spin==-1;
    if any(ip), scatter(hitX{q}(ip)*1e3,hitZ{q}(ip)*1e3,3,'r','filled'); end
    if any(im), scatter(hitX{q}(im)*1e3,hitZ{q}(im)*1e3,3,'b','filled'); end
    xlabel('X, мм'); ylabel('Z, мм');
    title(sprintf('Детектор %d, Y=%.2f м',q,detectorY(q)));
    xlim(detectorHalfX*1e3*[-1 1]); ylim(detectorHalfZ*1e3*[-1 1]);
end
exportgraphics(fig,fullfile(outDir,'07_detector_scatter_3screens.png'),'Resolution',220);

%% 11. ГИСТОГРАММЫ ПО Z НА ТРЁХ ДЕТЕКТОРАХ
fig=figure('Color','w','Name','Гистограммы детекторов');
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
edges=linspace(-detectorHalfZ,detectorHalfZ,121)*1e3;
for q=1:3
    nexttile; hold on; grid on;
    valid=isfinite(hitZ{q});
    zp=hitZ{q}(valid & spin==1)*1e3;
    zm=hitZ{q}(valid & spin==-1)*1e3;
    if ~isempty(zp), histogram(zp,edges,'Normalization','pdf','FaceColor',[0.85 0.2 0.2],'FaceAlpha',0.45,'EdgeColor','none'); end
    if ~isempty(zm), histogram(zm,edges,'Normalization','pdf','FaceColor',[0.2 0.35 0.85],'FaceAlpha',0.45,'EdgeColor','none'); end
    xlabel('Z, мм'); ylabel('Плотность, 1/мм');
    title(sprintf('Детектор %d',q));
    legend('\mu=+\mu_B','\mu=-\mu_B','Location','best');
end
exportgraphics(fig,fullfile(outDir,'08_detector_histograms.png'),'Resolution',220);

%% 12. ВЛИЯНИЕ ТОКА НА ОСЕВОЕ ПОЛЕ (раздел 5.2)
fig=figure('Color','w','Name','Влияние тока');
hold on; grid on;
currentFactors=[0.5 1 1.5];
for q=1:3
    Iq=p.I*currentFactors(q);
    plot(yGrid,Baxis*(Iq/p.I),'LineWidth',1.5,'DisplayName',sprintf('I = %.0f А',Iq));
end
xlabel('Y, м'); ylabel('B_y(0,Y), Тл');
title('Влияние силы тока на осевую компоненту поля');
legend('Location','best');
exportgraphics(fig,fullfile(outDir,'09_axis_field_different_currents.png'),'Resolution',220);

%% 13. СХЕМАТИЧЕСКОЕ СРАВНЕНИЕ С ОРИГИНАЛЬНЫМ ОПЫТОМ (раздел 5.3)
fig=figure('Color','w','Name','Сравнение с опытом Штерна—Герлаха');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; hold on; axis equal; grid on;
plot([0 0],[-1 1],'k','LineWidth',2);
plot([-0.15 0.15],[-0.45 0.45],'k','LineWidth',2);
plot([-0.15 0.15],[0.45 -0.45],'k','LineWidth',2);
text(0.05,0.82,'Две компоненты');
title('Схема характерного результата опыта');
xlabel('Условная координата'); ylabel('Поперечное отклонение');
nexttile; hold on; grid on;
for q=1:3
    valid=isfinite(hitZ{q});
    histogram(hitZ{q}(valid)*1e3,linspace(-80,80,100), ...
        'Normalization','count','DisplayStyle','stairs','LineWidth',1.2, ...
        'DisplayName',sprintf('Детектор %d',q));
end
xlabel('Z, мм'); ylabel('Число атомов');
title('Распределения, рассчитанные в MATLAB'); legend('Location','best');
exportgraphics(fig,fullfile(outDir,'10_classic_vs_matlab.png'),'Resolution',220);

%% 14. СТАТИСТИКА, ТАБЛИЦЫ И СОХРАНЕНИЕ
Detector=(1:3)';
Y_m=detectorY(:);
Hits=zeros(3,1); MeanZ_plus_mm=nan(3,1); MeanZ_minus_mm=nan(3,1);
StdZ_plus_mm=nan(3,1); StdZ_minus_mm=nan(3,1);
for q=1:3
    valid=isfinite(hitZ{q});
    Hits(q)=nnz(valid);
    zp=hitZ{q}(valid & spin==1)*1e3;
    zm=hitZ{q}(valid & spin==-1)*1e3;
    if ~isempty(zp), MeanZ_plus_mm(q)=mean(zp); StdZ_plus_mm(q)=std(zp); end
    if ~isempty(zm), MeanZ_minus_mm(q)=mean(zm); StdZ_minus_mm(q)=std(zm); end
end
stats=table(Detector,Y_m,Hits,MeanZ_plus_mm,MeanZ_minus_mm,StdZ_plus_mm,StdZ_minus_mm);
writetable(stats,fullfile(outDir,'detector_statistics.csv'));

% Дополнительная таблица поля вдоль оси.
writematrix([yGrid(:),Baxis(:),dBaxis(:)],fullfile(outDir,'axis_field_and_gradient.csv'));

results=struct();
results.parameters=p;
results.N_atoms=N_atoms; results.temperature_K=T;
results.vMeanMaxwell_m_s=vMeanMB; results.speedRange_m_s=[vMin vMax];
results.initialSpeed_m_s=speed0; results.spin=spin;
results.slit1=struct('Y_m',slit1Y,'radius_m',slit1Radius);
results.slit2=struct('Y_m',slit2Y,'halfWidth_m',slit2HalfWidth, ...
    'halfHeight_m',slit2HalfHeight,'zShift_m',slit2Zshift, ...
    'shape','intersection of two parabolas');
results.rGrid_m=rGrid; results.yGrid_m=yGrid;
results.By_T=By; results.Br_T=Br; results.Bmag_T=Bmag;
results.dBydY_T_per_m=dBydY;
results.dBdr_T_per_m=dBdr; results.dBdy_T_per_m=dBdy;
results.detectorY_m=detectorY; results.hitX_m=hitX; results.hitZ_m=hitZ;
results.passedSlit1=nnz(passed1); results.passedSlit2=nnz(passed2);
results.statistics=stats;
save(fullfile(outDir,'stern_gerlach_results.mat'),'results','-v7.3');

fprintf('\nРасчёт завершён. Результаты: %s\n',outDir);
fprintf('Прошли первую щель: %d из %d\n',nnz(passed1),N_atoms);
fprintf('Прошли вторую щель: %d из %d\n',nnz(passed2),N_atoms);
fprintf('Столкновения и рассеяние на остаточном газе не моделируются.\n');
disp(stats);
fprintf(['Примечание: искусственное разделение не добавляется. ' ...
    'Реальная величина разделения определяется рассчитанным полем и параметрами модели.\n']);
end
