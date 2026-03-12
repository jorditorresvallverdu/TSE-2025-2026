#Author: Jordi Torres Vallverdú
#Start Date: 12/03/2025

###Main

#Nunmber of firms and states 
const max_firms= 3;
const kmax= 19;
const start_firms= 1;

#entry and exit of firms 
const entry_low= 0.15;
const entry_high= 0.25;
const scrap_value= 0.1;
const entry_at= 4; 
const beta = 0,925; 
const delta= 0.7;

#investment


% entry-exit
c.ENTRY_LOW = 0.15; % uniform distribution support: lower bound
c.ENTRY_HIGH = 0.25; % uniform distribution support: higher bound
c.SCRAP_VAL = 0.1; % scrap value
c.ENTRY_AT = 4; % efficiency level at which new firms enter
c.BETA = 0.925; % discount factor
c.DELTA = 0.7; % probability of industry aggregate decline

% investment cost
c.INV_MULT = 3; % investment cost parameter

% profit
c.INTERCEPT = 3; % cournot demand intercept 
c.FIXED_COST = 0.2; % cournot fixed cost
c.GAMMA = 1; % cournot marginal cost coefficient

% convergence tolerence:
c.TOL = 0.1; % tolerance

% indicators and name prefixes that make it easier to interpret saved
% results
c.PROFIT_DONE = 0; % indicator for having finished computing profit
c.EQL_DONE = 0; % indicator for having finihsed equilibrium computation
c.PREFIX = 'cc'; % prefix representing Cournot competition in saved results

% simulation
c.DS_WSTART = [c.ENTRY_AT+2; zeros(c.MAX_FIRMS-1,1)]; % initial state for simulation
c.DS_NSIMX = 10000; % number of simulation periods

%% Test decode function
%%% Task %%%%%

%%%%%%%%%%%%%%

%% Baseline
% Compute static profit:
static_profit(c);

% Solve dynamic equilibrium:
eql_ma(c);

% Simulate entry & exit:
ds_ma(c,'baseline.tex');

%% Low entry cost
% new entry cost parameters
c.ENTRY_LOW = 0.01;
c.ENTRY_HIGH = 0.11;

% Compute static profit:
static_profit(c);

% Solve dynamic equilibrium:
eql_ma(c);

% Simulate entry & exit:
ds_ma(c,'low_entry_cost.tex');

