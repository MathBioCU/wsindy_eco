%% load data
dr = '~/Dropbox/Boulder/research/data/dukic collab/';
% df = 'trap_data10.csv';
df = 'trap_idaho_1977-2023.csv';
dat = readtable([dr,df]);% [a,~]=unique(dat.PlotName);
num_years = 47;% 39
x_nan = reshape(dat.trap_mean,num_years,[]);
t = 0:num_years-1;
latlon = [dat.lat(1:num_years:end) dat.lat(1:num_years:end)];
[aa,bb]= ndgrid((1:length(latlon)),(1:length(latlon)));
dd = arrayfun(@(i,j) distance(latlon(i,:),latlon(j,:)),aa(:),bb(:));
dd = reshape(dd,length(latlon),[]);
x = x_nan;
df = @(d)d.^2;
for j=1:size(x,1)
    for i=1:size(x,2)
        if isnan(x(j,i))
            v = x(j,:);
            d = 1./df(max(dd(i,:),0.01));
            d = d(~isnan(v));
            v = v(~isnan(v));
            d = d/sum(d);
            if any(isinf(d))
                disp([v' d'])
            end
            x(j,i) = dot(v,d);
        end
    end
end

i=1;
plot([x(:,i) x_nan(:,i)],'o-')

%% weather data
df = 'ym_weather_idaho_1979-2023.csv';
t_dat = readtable([dr,df]);% [a,~]=unique(dat.PlotName);
num_years = 47;% 39
keys = unique(t_dat.key);
t_x = cell(1,length(keys));
for j=1:length(keys)
    t_x{j} = t_dat.sum_pr(ismember(t_dat.key,keys(1)));
end
t_x = movmean(cell2mat(t_x),12);
t_x = t_x(1:12:end,:);
t_x = [t_x(1:2,:);t_x];

%% examine corr
ind = 8;
plot([t_x(:,ind)*400 x_nan(:,ind)])

%% interpolate onto finer grid 

Nq = num_years;
tq = linspace(t(1),t(end),Nq);
y = interp1(t,x,tq,'cubic');
y = max(y,0);
t = tq;
M = length(t);

%% cluster trajectories
N = 1; % num clusters
a = kmeans(y',N);
[c,d]=max(arrayfun(@(i)sum(a==i),1:N));
x = x(:,a==d);
ntraj = size(x,2);

%% view data
plot(t,y)

%% get data object
Uobj = arrayfun(@(i)wsindy_data(x(:,i),t(:)),1:ntraj);

%% add noise
noise_ratio = 0.; rng(1); rng_seed = rng().Seed; rng(rng_seed);
Uobj.addnoise(noise_ratio,'seed',rng_seed);

%% get lib tags
nstates = Uobj.nstates;
polys = [0:2];
trigs = [];
tags_1 = get_tags(polys,[],nstates);
tags_2 = get_tags([],trigs,nstates);
tags = prodtags(tags_1,tags_2);
% tags = [tags;{@(x)x./(0.1+x.^2);@(x)x./(0.01+x.^2);@(x)x./(0.001+x.^2)}];
lib = library('tags',tags);
% lib.add_terms({term('fHandle',@(x)x.^2,'linOp',1),term('fHandle',@(x)x,'linOp',1)});

%% get test function

phifuns = {optTFcos(2,0),optTFcos(2,1),optTFcos(2,2)};
tf = {};
for i=1:length(phifuns)
    tf = [tf;arrayfun(@(U){testfcn(U,'stateind',1,'phifuns',phifuns(i),'meth','FFT','param',2)},Uobj','uni',0)];
end
WS = wsindy_model(repmat(Uobj(:),length(phifuns),1),lib,tf,'lhsterms',{term('ftag',1,'linOp',2)});

optm = WS_opt();
[WS,loss_wsindy,its,G,b] = optm.MSTLS(WS,'alpha',0.01);
% [WS,w_its,res,res_0,CovW] = optm.wendy(WS,'diag_reg',10^-6,'maxits',100);

%% 1,3,6,7  2,4,5 are irregular
for traj_ind = 1
if traj_ind>0
    [PKS,LOCS] = findpeaks(Uobj(traj_ind).Uobs{1});
    [a,b] = sort(PKS,'descend');

    w_plot = WS.weights;
    d = max(PKS);    % d = mean(a(a>10));
    c = -d^2*w_plot(cellfun(@(tt)and(isequal(tt.ftag,2),isempty(tt.linOp)),lib.terms))/3 - ...
        d*w_plot(cellfun(@(tt)and(isequal(tt.ftag,1),isempty(tt.linOp)),lib.terms))/2;
    w_plot(cellfun(@(tt)and(isequal(tt.ftag,0),isempty(tt.linOp)),lib.terms)) = c;
    rhs_learned = WS.get_rhs('w',w_plot);

    tol_NLS = 10^-6;
    options_NLS = odeset('RelTol',tol_NLS,'AbsTol',tol_NLS*ones(1,nstates));
    t_learned = t;
    x0_reduced = [Uobj(traj_ind).Uobs{1}(1);diff(Uobj(traj_ind).Uobs{1}(1:2))];
    tp = LOCS(b(1:3));
    tp = LOCS(b);
    tp = ':';
    % fun = @(x0) norm(sol(x0,rhs_learned,t_learned,options_NLS,t(tp))-Uobj(traj_ind).Uobs{1}(tp));
    % x0_reduced = fmincon(fun,x0_reduced,[],[],[],[],[0 0]);
    fun = @(x) norm(sol([x 0],rhs_learned,t_learned,options_NLS,t(tp))-Uobj(traj_ind).Uobs{1}(tp));
    x0_reduced = [fmincon(fun,x0_reduced(1),[],[],[],[],0,d) 0];
    % x0_reduced = [0 0 ];

    tol_dd = 10^-12;
    options_ode_sim = odeset('RelTol',tol_dd,'AbsTol',tol_dd*ones(1,nstates));
    t_learned = linspace(t(1),t(end),2000);
    [t_learned,xH0_learned]=ode15s(@(t,x)rhs_learned(x),t_learned,x0_reduced,options_ode_sim);

    figure(7);clf
    Uobj(traj_ind).plotDyn;hold on
    plot(t_learned,xH0_learned(:,1),'--')
    hold off
    title(x0_reduced)
    drawnow
    
    % figure(8)
    % if exist('w_its','var')
    %     plot_wendy;
    % end
end
end

function x = sol(x0,rhs,t,opt,tq)
    [t_l,x]=ode15s(@(t,x)rhs(x),t,x0,opt);
    try
        x = interp1(t_l,x(:,1),tq);
        x = fillmissing(x,'constant',0);
    catch
        x = repmat(x0(:,1),length(tq),1);
    end
end