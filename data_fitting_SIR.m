% Do some housekeeping
clc
clear
close all

% Data originating from Johns Hopkins University
fileName = [tempdir 'time_series_covid19_confirmed_global.csv'];
if ~isfile(fileName)
    url = "https://github.com/CSSEGISandData/COVID-19/raw/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv";
    fileName = websave(fileName, url);
end
opts = detectImportOptions(fileName); % Detect import parameters

% Fix range and delimiter
opts.DataLines = [2, 5];
opts.Delimiter = ",";

% Fix column names and types for first columns
opts.VariableNames(1:4) = [{'ProvinceState'}, {'CountryRegion'}, {'Lat'}, {'Long'}];
opts.VariableTypes(1:4) = [{'string'}, {'categorical'}, {'double'}, {'double'}];

% Fix file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Fix variable properties
opts = setvaropts(opts, "ProvinceState", "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["ProvinceState", "CountryRegion"], "EmptyFieldRule", "auto");

% Import the data
data = readtable(fileName, opts);

%% calculations
avg = 7; % Days to average
c = size(data, 2); % Number of columns
ds = c-4; % Number of days with data
t0 = datetime('22/1/2020'); % Date of first data
lastDate = t0 + days(ds-1); % Date of last data

% Prepare variables for cumulate country data
cats = unique(data.CountryRegion); % Countries
lc = length(cats); % Number of countries
data2 = data; % Copy data to force same structure for both dataset
data2(lc+1:end, :) = []; % Remove unnecessary rows

% Cumulate country data
for i = 1:lc
    tmp1 = mean(data{data.CountryRegion == cats(i), 3:4}, 1); % Coordinates - mean
    tmp2 = sum(data{data.CountryRegion == cats(i), 5:end}, 1); % Cases - sum
    data2(i, 1:2) = table("", cats(i)); % Assign country to first columns
    data2(i, 3:end) = array2table([tmp1 tmp2]); % Assign mean of coordiantes and sum of cases
end
%data2
%clear tmp1 tmp2
%% smoothing data
cc = c+ds-avg+1;
data2{:, c+1:cc} = movmean(data2{:, 5:c}, avg, 2, 'Endpoints', 'discard');

% Calculate daily differences
data2{:, cc+1:cc+ds-avg} = diff(data2{:, c+1:cc}')';
%size(data2)
%%
% Filter by region
% catsTerritory = {'Bosnia and Herzegovina', 'Bulgaria', 'Croatia', 'Macedonia', 'Hungary', 'Montenegro', 'Romania', 'Serbia', 'Slovenia'};
% catsEurope = {'Italy', 'Spain', 'Germany', 'France', 'UK', 'Switzerland', 'Belgium', 'Netherlands', 'Austria'};
% catsG7 = {'Canada','France', 'Germany', 'Italy', 'Japan', 'United Kingdom', 'US'};
% filter = catsTerritory;
% data3 = data2(ismember(data2.CountryRegion, filter), :); % Data of selected countries
% X = data3{:, c+2:cc}; % Number of Total Confirmed cases
% Y = data3{:, cc+1:end}; % New confirmed cases
% L = cellstr(data3.CountryRegion); % Country list for legend
% [ctrs, dys] = size(X); % Helper variables

X1=data2{:, c+2:cc};% Number of Total Confirmed cases
Y1= data2{:, cc+1:end};% New confirmed cases
L1 = cellstr(data2.CountryRegion);
[ctrs1, dys1] = size(X1);

loglog(X1',Y1'); % Create plot with logarithmic scale
formatOut = 'mm/dd';
lastDate = datestr(lastDate, formatOut);
title(['Trajectory of COVID-19 Confirmed Cases until ', char(lastDate)]);
xlabel('Total Confirmed Cases');
ylabel('New Confirmed Cases (in the past week)');
legend(L1, 'Location', 'best');
%% fitting the data to the model to find the nearest beta and gamma
% 
%for day 1 the initial number of infected was 0 and this is an equilibrium
% state for our SIR model. So we chose a different initial day with positive 
% value for infected number of people
t0 = 67;

%t_end denotes the final time point for the model
t_end = cc-c-t0+1;

%population of each country taken from Worldometer
population = table([38928346;2877797;43851044;77265]);

%total number of countries
total = 4;

%fitting using fminsearch and the SIR model with constant parameters
for i = 1:total

    %data is taken from t0 to the last available day
    inf_data = data2{i, c+t0:cc};

    %initial condition is taken as the data at t0
    InC =data2{i, c+t0};
    N = population.Var1(i); 
    obj = @(params) obj_func(params(1),params(2),N,inf_data,t_end,InC);%params1 = beta params2 = gamma
    [params_fit,feval] = fminsearch(obj,[0.3,0.1]');
    beta(i) = params_fit(1);
    gamma(i) = params_fit(2);
    Inf = SIR_model_1(beta(i),gamma(i),N,t_end,InC);
    figure
        plot(t0:cc-c,Inf,'LineWidth',3)
        hold on
        plot(t0:cc-c,inf_data,'LineWidth',3)
        legend('SIR Model','data')
        title(data2.CountryRegion(i))
        set(gca,'FontSize',12,'FontWeight','bold')
        xlabel('days')
        ylabel('number of infected people')
end
%% Improving fit
% In this section we are using a time dependent function for beta and this
% is directly involving the data and the daily differences. We use this
% time-dependent beta in the model which inturn is dependent on the
% recovery rate gamma. We write SIR_time_model as the function for the
% differential equations and the objective function is the quadratic error
% function. 
for i = 1:total
inf_data_1 = data2{i, c+t0:cc};
daily_diff = data2{i, cc+t0:cc+ds-avg};
InC1 =data2{i, c+t0};
N1 = population.Var1(i);
%params1 = beta params2 = gamma
obj = @(params) obj_func_new(params,N1,inf_data_1,daily_diff,t_end,InC1);
[params_fit,feval] = fminsearch(obj,0.1);
Inf = SIR_time_model(inf_data_1,daily_diff,params_fit,N1,t_end,InC1);
    figure
        plot(1:length(Inf),Inf,'o','LineWidth',1)
        hold on
        plot(1:length(Inf),inf_data_1(2:end-1),'LineWidth',3)
        legend('SIRModel','data')
        title(data2.CountryRegion(i))
        set(gca,'FontSize',12,'FontWeight','bold')
        xlabel('days')
        ylabel('number of infected people')

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%%Objective function

function cost = obj_func(beta,gamma,N,inf_data,tend,init_I)


Inf = SIR_model_1(beta,gamma,N,tend,init_I);

%cost_functional is the square of the difference of the solution from the model and
%the real data
cost = norm(Inf - inf_data)^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function for the evolution of disease using SIR model with constant
%parameters
function Inf = SIR_model_1(beta,gamma,N,tend,init_I)
% X be the unknown variable that is evolving with time. X = [S;I;R]. We
% will use ode45 solver to find the solution for the system given by 
% X' = [-beta*I*S/N;beta*I*S/N -gamma*I; gamma*I];
% We need to use a total population N
dt = 1;
t0 = 1;
%gamma = 1/10;
InC = [N-init_I,init_I,0]';
func = @(t,X) [-beta*X(2)*X(1)/N;beta*X(2)*X(1)/N-gamma*X(2); gamma*X(2)];
[t,y] = ode45(func,t0:dt:tend,InC);

% time span length
l = length(t);

% The infected numbers
Inf = y (l+1:2*l);

% figure
% plot(t,y)
% legend('S','I','R')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function for the upgraded SIR model with a time-dependent beta as the
%parameter and gamma is to be fitted using fminsearch
function Inf = SIR_time_model(data,diff,gamma,N,tend,init_I)

dt = 1;
t0 = 1;
%gamma = 1/10;
% taking only upto tend-1 so that we have same length of data and data_diff
tspan = t0:dt:tend-1; 
% time span length
l = length(tspan);
beta_t = @(k,S,gamma) (N/S)*((1/data(k))*(diff(k)/dt)-gamma);

InC = [N-init_I,init_I,0]';
func = @(t,X) [-beta_t(t,X(1),gamma)*X(2)*X(1)/N;beta_t(t,X(1),gamma)*X(2)*X(1)/N-gamma*X(2); gamma*X(2)];
%[~,y] = ode45(func,tspan,InC);

% % The infected numbers
% Inf = y (l+1:2*l);

%explicit euler solver
y(:,1) = InC;
for n = 2:l
    y(:,n) = y(:,n-1)+dt*func(tspan(n-1),y(:,n-1));
end
% The infected numbers
Inf = y (2,2:l);
end


% the new cost_function
function cost = obj_func_new(gamma,N,inf_data,diff,tend,init_I)
Inf = SIR_time_model(inf_data,diff,gamma,N,tend,init_I);
%cost_functional is the square of the difference of the data from the model and
%the real data
cost = norm(Inf - inf_data(2:end-1))^2;
end





