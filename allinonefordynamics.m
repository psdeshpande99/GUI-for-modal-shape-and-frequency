%% Data Input and Processing
prompt1 = {'Enter no of Degrees of Freedom'};
dfin  = {'1'};
n = inputdlg(prompt1,'Input',[1,60],dfin);
n = str2double(n);
response = zeros (n,n);
response2 = zeros (n,n);
M = zeros(n,n);
for i=1:n
    c = i:n;
    c = arrayfun(@num2str, c, 'uni',0);
    prompt2 = c;
    inbuiltip = zeros(1,n-i+1);
    inbuiltip = arrayfun(@num2str, inbuiltip, 'uni',0);
    dfin2  = inbuiltip;
    dlgtitle = (['for mass ',num2str(i),' th']);
    mat = inputdlg(prompt2,dlgtitle,[1,60],dfin2);
    response(i,i:n) = str2double(mat');
end

c = 1:n;
c = arrayfun(@num2str, c, 'uni',0);
prompt2 = c;
inbuiltip = zeros(1,n);
inbuiltip = arrayfun(@num2str, inbuiltip, 'uni',0);
dfin2  = inbuiltip;
dlgtitle = ('for mass matrix');
mat = inputdlg(prompt2,dlgtitle,[1,60],dfin2);
response2(1,:) = str2double(mat');

c = 1:2;
c = arrayfun(@num2str, c, 'uni',0);
prompt3 = c;
inbuiltip = zeros(1,2);
inbuiltip = arrayfun(@num2str, inbuiltip, 'uni',0);
dfin3  = inbuiltip;
dlgtitle = ('Multipliers');
mat = inputdlg(prompt3,dlgtitle,[1,60],dfin3);
response5(1,:) = str2double(mat');


c = 1:2;
c = arrayfun(@num2str, c, 'uni',0);
prompt3 = c;
inbuiltip = zeros(1,2);
inbuiltip = arrayfun(@num2str, inbuiltip, 'uni',0);
dfin3  = inbuiltip;
dlgtitle = ('ws required');
mat = inputdlg(prompt3,dlgtitle,[1,60],dfin3);
response4(1,:) = str2double(mat');

response = response5(1,1).*response;
response2 = response5(1,2).*response2;
for j=1:n
    M(j,j) = response2(1,j);
end

K = -1.*response;

for j=1:n
    for k=1:n
        if(K(k,j)~=0)
        K(j,k) = K(k,j);
        end
    end
end
for i=1:n
    K(i,i) = abs(sum(K(i,:),2));
end


%% Eigen solver method
[modeeigen,w_squareeigen]=eig(K,M);
freqeigen=sqrt(diag(w_squareeigen));
for i=1:n
    modeeigen(:,i) = modeeigen(:,i)./modeeigen(1,i); 
end
    
    


%% Dunkerley Method %%
K_dunkerley = [response(1,1); diag(response,1)];
K_dunkerley = 1./K_dunkerley;
M_dunkerley = diag(M);
freqdunker = sqrt(1/sum((M_dunkerley.*cumsum(K_dunkerley))));

%% Stodala Method
K_stodala = [response(1,1); diag(response,1)];
ad=ones(n,1);% Asssumed deflection
for u=1:100000
    fi = ad.*M_dunkerley;
    fs = flip(cumsum(flip(fi)));
    sd = cumsum(fs./K_stodala);
    modestodala = sd./sd(1);
    freqstodala = sqrt(sum(ad)./sum(sd));
    ad = modestodala;
    
end

%% Influence Method
x=ones(n,1);
for i=1:1000
    z = inv(K)*M*x;
    r = z(1);
    modeinfluence = z./z(1);
    if(modeinfluence == x)
        freqinfluence(1,1) = sqrt(1/r);
        break
    else
        x = modeinfluence;
    end
end
modeinfluence(:,n) = modeinfluence;
modeinfluence(:,2:n) = 1;
massmultiplier = modeinfluence .* M_dunkerley;
for i=1:n-1
    multiplier = zeros(n,n);
    for j=1:n
        multiplier(j,j) = 1;
    end
    multiplier(i,:) = -1*(massmultiplier(:,i)./massmultiplier(1,1));
    multiplier(1,1) = 0;
    for p=1:1000
    z = inv(K)*M*multiplier*x;
    r = z(1);
    modeinfluence(:,i+1) = z./z(1);
    if(modeinfluence(:,i+1) == x)
        freqinfluence(1,i+1) = sqrt(1/r);
        break
    else
        x = modeinfluence(:,i+1);
    end
    freqinfluence(1,i+1) = sqrt(1/r);
    end
    massmultiplier = modeinfluence .* M_dunkerley;
end
    

%% Damping Matrix Rayleigh
newfreqeigen = freqeigen';
zeta = zeros(1,n);
zeta(1,response4(1,1)) = 0.015;
zeta(1,response4(1,2)) = 0.015;
for i=1:10000
    a1 = 2*zeta(1,response4(1,1))/(newfreqeigen(response4(1,1)) + newfreqeigen(response4(1,2)));
    a0 = 2*zeta(1,response4(1,1))*(newfreqeigen(response4(1,1)) * newfreqeigen(response4(1,2)))/(newfreqeigen(response4(1,1)) + newfreqeigen(response4(1,2)));
    for j=1:n
        if(response4(1,1)==j || response4(1,2)==j)
            continue
        else
            zeta(1,j) = (a0/(2*newfreqeigen(1,j))) + (newfreqeigen(1,j)*a1)/2;
        end
    end
    
end
C = a0*M + a1*K;


%% Damping by Caughey
newfreqeigen = freqeigen';
zetacaughey = zeros(n,1);
zetacaughey(1:n,1) = 0.015;
Multiplycaughey = zeros(n,n);
Ccaughey = zeros(n,n);
for i=1:n
    for j=0:n-1
        Multiplycaughey(i,j+1) = (newfreqeigen(i)^(2*j-1))/2;
    end
end
X = linsolve(Multiplycaughey,zetacaughey);
for i=1:n
    Ccaughey = M*X(i)*(inv(M)*K)^(i-1) + Ccaughey;
end




%% Final display
%uiwait(msgbox([ ' The Matrix is ', num2str(K) ]));
response_update = diag(response);
response_new = response_update(2:n);

n_array = 1:n;
f = figure('Name','Lumped mass all method compare','NumberTitle','off','Units','normalized','Position',[0 0 0.4 0.4]);
uit = uitable('Data', K,'BackgroundColor',[1 0 1],'FontSize',10);
uit.Position = [20 320 450 250];
set(uit,'ColumnWidth',{75})
txt_title1 = uicontrol('Style', 'text', 'Position', [50 570 200 20], 'String', 'Stiffness Matrix (N/m)','ForegroundColor',[0.8500 0.3250 0.0980],'FontSize',12);
uit = uitable('Data', M,'BackgroundColor',[1 0 1],'FontSize',10);
uit.Position = [20 10 450 250];
set(uit,'ColumnWidth',{75})
txt_title2 = uicontrol('Style', 'text', 'Position', [50 270 200 20], 'String', 'Mass Matrix (kg)','ForegroundColor',[0.8500 0.3250 0.0980],'FontSize',12);
uit = uitable('Data', freqeigen','BackgroundColor',[0.9290 0.6940 0.1250],'FontSize',10);
uit.Position = [480 530 350 38.5];
set(uit,'ColumnWidth',{75})
txt_title3 = uicontrol('Style', 'text', 'Position', [450 570 400 20], 'String', 'All Frequencies By Eigen Solver (rad/s)','ForegroundColor','r','FontSize',12);
uit = uitable('Data', modeeigen,'BackgroundColor',[0.9290 0.6940 0.1250],'FontSize',10);
uit.Position = [480 320 350 180];
set(uit,'ColumnWidth',{75})
txt_title4 = uicontrol('Style', 'text', 'Position', [470 500 300 20], 'String', 'Mode shapes by Eigen Solver','ForegroundColor','r','FontSize',12);
uit = uitable('Data', freqdunker','BackgroundColor',[0.9290 0.6940 0.1250],'FontSize',10);
uit.Position = [480 220 120 38.5];
set(uit,'ColumnWidth',{75})
txt_title5 = uicontrol('Style', 'text', 'Position', [470 270 150 20], 'String', 'Wn by Dunkerley','ForegroundColor','r','FontSize',12);
uit = uitable('Data', freqstodala','BackgroundColor',[1 1 0],'FontSize',10);
uit.Position = [480 150 120 38.5];
set(uit,'ColumnWidth',{75})
txt_title6 = uicontrol('Style', 'text', 'Position', [470 190 150 20], 'String', 'Wn by Stodala','ForegroundColor',[0 0 0],'FontSize',12);
uit = uitable('Data', modestodala,'BackgroundColor',[1 1 0],'FontSize',10);
uit.Position = [660 10 100 180];
set(uit,'ColumnWidth',{75})
txt_title7 = uicontrol('Style', 'text', 'Position', [640 190 150 20], 'String', 'Mode by Stodala','ForegroundColor',[0 0 0],'FontSize',12);
uit = uitable('Data', zetacaughey,'BackgroundColor',[0.3010 0.7450 0.9330],'FontSize',10);
uit.Position = [850 530 350 38.5];
set(uit,'ColumnWidth',{75})
txt_title8 = uicontrol('Style', 'text', 'Position', [790 570 500 20], 'String', 'Damping Ratios of all modes by Caughey','ForegroundColor','b','FontSize',12);
uit = uitable('Data', Ccaughey,'BackgroundColor',[0.3010 0.7450 0.9330],'FontSize',10);
set(uit,'ColumnWidth',{75})
uit.Position = [850 320 350 180];
txt_title9 = uicontrol('Style', 'text', 'Position', [830 500 400 20], 'String', 'Damping Matrix by Caughey Method','ForegroundColor','b','FontSize',12);
subplot(2,3,6)
for i=1:n
    plot(modeeigen(:,i),1:n,'x-','DisplayName',num2str(n_array(i)))
    hold on;
end
title('Mode shape plot')
grid on;
ylim([1 n]);
xlabel('displacements');
ylabel('mass number');
legend show





