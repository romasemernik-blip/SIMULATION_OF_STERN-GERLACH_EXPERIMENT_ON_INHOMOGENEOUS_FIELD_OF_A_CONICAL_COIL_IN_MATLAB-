function field_conical_coil_report_method()
% Поле конической катушки по алгоритму, описанному в разделе 3.2 реферата.
% Ось Y совпадает с осью катушки и направлением атомного пучка.
% Метод: сумма осевых полей кольцевых слоёв; центральная разность dBy/dY;
% приближение Br = -(r/2)*dBy/dY; затем |B| = hypot(Br,By).
% Важно: это осесимметричное приближение вблизи оси, а не точное поле
% реальной винтовой обмотки вблизи проводников.

clc; close all;
mu0 = 4*pi*1e-7;       % Гн/м
n = 3200;              % плотность витков, 1/м
I = 650;               % ток, А
a = 0.32;              % входной радиус, м
b = 0.44;              % выходной радиус, м
L = 1.15;              % длина, м
y0 = 0.0;              % начало катушки, м
Nlayers = 720;         % узлы численного суммирования

y = linspace(-0.35,2.30,650);
r = linspace(0,0.55,300);
s = linspace(y0,y0+L,Nlayers);
ds = s(2)-s(1);
w = ds*ones(size(s)); w([1 end]) = ds/2;
R = a + (b-a)*(s-y0)/L;

[RR,YY] = ndgrid(r,y);
By = zeros(size(RR));
for k = 1:Nlayers
    % Осевое поле одного кругового витка, умноженное на n*ds.
    dBy = mu0*n*I*w(k)*R(k)^2 ./ ...
        (2*(R(k)^2 + (YY-s(k)).^2).^(3/2));
    By = By + dBy;
end

dy = y(2)-y(1);
dBydY = zeros(size(By));
dBydY(:,2:end-1) = (By(:,3:end)-By(:,1:end-2))/(2*dy);
dBydY(:,1) = (By(:,2)-By(:,1))/dy;
dBydY(:,end) = (By(:,end)-By(:,end-1))/dy;

Br = -0.5 .* RR .* dBydY;
Bmag = hypot(Br,By);
dByAxis = gradient(By(1,:),dy);

figure('Color','w','Position',[100 100 1300 800]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile;
contourf(y,r,Bmag,36,'LineColor','none'); colorbar;
hold on; plot([y0 y0+L],[a b],'w-','LineWidth',2);
plot([y0 y0+L],[-a -b],'w-','LineWidth',2);
xlabel('Y, м'); ylabel('r, м');
title('Модуль магнитной индукции |B|');

nexttile;
contourf(y,r,By,36,'LineColor','none'); colorbar;
hold on; plot([y0 y0+L],[a b],'w-','LineWidth',2);
plot([y0 y0+L],[-a -b],'w-','LineWidth',2);
xlabel('Y, м'); ylabel('r, м');
title('Продольная компонента B_y');

nexttile;
plot(y,By(1,:),'LineWidth',1.8); grid on;
xlabel('Y, м'); ylabel('B_y(0,Y), Тл');
title('Осевое поле вдоль оси катушки');

nexttile;
plot(y,dByAxis,'LineWidth',1.8); yline(0,'k-'); grid on;
xlabel('Y, м'); ylabel('dB_y/dY, Тл/м');
title('Градиент продольной компоненты на оси');

sgtitle('Магнитное поле конической катушки — методика раздела 3.2');
exportgraphics(gcf,'Поле_конической_катушки_по_реферату.png','Resolution',220);

save('Поле_конической_катушки_по_реферату.mat', ...
    'y','r','By','Br','Bmag','dByAxis','mu0','n','I','a','b','L','y0');
end
