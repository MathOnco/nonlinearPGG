function [P,F,fraction,t,Pvec,Fvec,groupsize,nonemptyCompartment] = multi_compartment_model(options)

close all;

Ncomp         = 70;                  % Number of compartments
t0            = 0;                    % Initial time        
tf            = 1e9;                  % Final time
Ts            = min(0.1*tf,100);       % Increment time (where we redistribute)
b             = 1;                    % Production/Consumption rate
delta_P       = 0.1;                  % Death rate of producers
delta_F       = 0.1;                  % Death rate of free-riders
P0            = 50;                  % Initial population size
F0            = 50;                   % Free-rider
K             = 1e3;                  % Carrying capacity
tol           = 1e-2;

plotting   = false;
final_plot = false;
random     = 'multinomial';
hysteresis = false;
exitcriterion = true;
shuffling = true;

mixing = false;      % If true, shuffle after Ts

% Growth function parameters
alpha = 1;
sigma = 3;
beta  = 2;
kappa = 0.5;

if ~exitcriterion
    tf = 3e3;
    if ~shuffling; Ts = 1e9; end
end

if nargin == 1
    if isfield(options,'Ncomp'); Ncomp = options.Ncomp; end
    if isfield(options,'alpha'); alpha = options.alpha; end
    if isfield(options,'beta'); beta = options.beta; end
    if isfield(options,'sigma'); sigma = options.sigma; end
    if isfield(options,'delta_P'); delta_P = options.delta_P; end
    if isfield(options,'delta_F'); delta_F = options.delta_F; end
    if isfield(options,'kappa'); kappa = options.kappa; end
    if isfield(options,'K'); K = options.K; end
    if isfield(options,'Ts'); Ts = options.Ts; end
    if isfield(options,'tf'); tf = options.tf; end
    if isfield(options,'P0'); P0 = options.P0; end
    if isfield(options,'F0'); F0 = options.F0; end
    if isfield(options,'plotting'); plotting = options.plotting; end
end

changevar = 'Ncomp';

% Growth functions
p1 = @(x,pro) (1 + (pro-1)*x)/pro;
p2 = @(x,pro) (pro-1)*x/pro;

G1 = @(x,pro,b) alpha*(1+exp(sigma))/(1+exp(sigma-b*p1(x,pro))) - kappa;
G2 = @(x,pro,b) alpha*(1+exp(sigma))/(1+exp(sigma-b*p2(x,pro)));

% Initialize
y = zeros(1,2*Ncomp);

switch random
    case 'multinomial'
    prob = ones(Ncomp,1)/Ncomp;
    y(1:2:end) = mnrnd(Ncomp,P0,prob);
    y(2:2:end) = mnrnd(Ncomp,F0,prob);
    case 'uniform'
        PickNewP = rand(Ncomp,1);
        DensityP = PickNewP/sum(PickNewP);
        PickNewF = rand(Ncomp,1);
        DensityF = PickNewF/sum(PickNewF);
        y(end,1:2:end) = P0*DensityP;
        y(end,2:2:end) = F0*DensityF;
    otherwise
        y(1:2:end) = P0/Ncomp;
        y(2:2:end) = F0/Ncomp;
end

t = t0;
breakout = false;
count = 0;
counter = 1;

medianPtotal = 1e10;
medianFtotal = 1e10;
GoDownOnly = false;
output = 1;
fraction = [];
thold = 0;

while t0 < tf
    
    if exist('degofchange1b','var')
        degofchange1a = degofchange1b;
        degofchange2a = degofchange2b;
    else
        degofchange1a = 1e10;
        degofchange2a = 1e10;
    end
    tspan = [t0 min(t0+Ts,tf)];
    [tnew,ynew,te] = ode45(@(t,y)RHS(t,y,Ncomp,b,delta_P,delta_F,K,G1,G2,beta),tspan,y(end,:));
    
    % Append to solution and time vector
    y = [y;ynew];
    t = [t;tnew];
    
    % Move new starting time forward
    t0 = t0+Ts;
    
    % Move around
    Ptotal(counter,1) = floor(sum(y(end,1:2:end)));
    Ftotal(counter,1) = floor(sum(y(end,2:2:end)));
    
    if t(end)-thold>1e3 && abs(medianPtotal/Ptotal(end)-1) < tol && ...
            abs(medianFtotal/Ftotal(end)-1) < tol
         count = count + 1;
%     elseif hysteresis && mod(counter,100) == 0
%         count = count - 1;
    end
    
    medianPtotal = median(Ptotal);
    medianFtotal = median(Ftotal);

    if exitcriterion && (Ptotal(counter)*Ftotal(counter) == 0 || count > 10)
        if hysteresis
            fprintf('%i\t%g\t%g\t%g\n',Ncomp,beta,Ptotal(counter)/K,Ftotal(counter)/K);
            fraction(output,:) = [beta Ptotal(counter)/(Ptotal(counter) + Ftotal(counter))];
            output = output + 1;
            if GoDownOnly || Ftotal(counter) == 0
                Ftotal(counter) = 20;
                switch changevar
                    case 'beta'
                        beta = beta-0.1;
                    case 'Ncomp'
                        Ncomp = Ncomp-1;
                        y = zeros(1,2*Ncomp);
                        prob = ones(Ncomp,1)/Ncomp;
                end
                GoDownOnly = true;
            else
                if Ptotal(counter)==0; Ptotal(counter) = 20; end
                switch changevar
                    case 'beta'
                        beta = beta+0.1;
                    case 'Ncomp'
                        Ncomp = Ncomp+1;
                        y = zeros(1,2*Ncomp);
                        prob = ones(Ncomp,1)/Ncomp;
                end
            end
            count = 0;
            Ptotal = Ptotal(counter,1);
            Ftotal = Ftotal(counter,1);
            counter = 1;
            tlast = t(end);
            if beta <= 0 || Ncomp < 1
                breakout = true;
            end
        else
            breakout = true;
            break;
        end
    end
    
    if breakout
        break;
    end

    if mixing
        switch random
            case 'multinomial'
                y(end,1:2:end) = mnrnd(Ncomp,Ptotal(counter),prob);
                y(end,2:2:end) = mnrnd(Ncomp,Ftotal(counter),prob);
            case 'uniform'
                PickNewP = rand(Ncomp,1);
                DensityP = PickNewP/sum(PickNewP);
                PickNewF = rand(Ncomp,1);
                DensityF = PickNewF/sum(PickNewF);
                y(end,1:2:end) = Ptotal*DensityP;
                y(end,2:2:end) = Ftotal*DensityF;
        end
    end
    

    if plotting && ~hysteresis
        drawnow;
        plot(t,[sum(y(:,1:2:end),2),sum(y(:,2:2:end),2)])
        xlabel('time'); ylabel('Population size');
        legend('Producers','Free-Riders','location','northwest');
        set(gca,'fontsize',16);
        pause(0.01);
    end
    
    counter = counter + 1;

end

groupsize = y(:,1:2:end) + y(:,2:2:end);

% Proportion of nonzero compartments (at least one individual >1)
nonemptyCompartment = sum(groupsize>1,2)/Ncomp;

if final_plot
    semilogy(t,[sum(y(:,1:2:end),2),sum(y(:,2:2:end),2)])
    xlabel('time'); ylabel('Population size');
    legend('Producers','Free-Riders','location','northwest');
    set(gca,'fontsize',16);
end

if breakout
    Pvec = floor(sum(y(:,1:2:end),2));
    Fvec = floor(sum(y(:,2:2:end),2));
    P    = Pvec(end);
    F    = Fvec(end);
else
    Pvec = sum(y(:,1:2:end),2);
    Fvec = sum(y(:,2:2:end),2);
    P = floor(mean(sum(y((t>0.9*t(end)),1:2:end),2)));
    F = floor(mean(sum(y((t>0.9*t(end)),2:2:end),2)));
end

% F = y(end,2);

end


function dydt = RHS(t,y,Ncomp,b,delta_P,delta_F,K,G1,G2,beta)

Y = sum(y);

for n = 1 : Ncomp

    P = y(2*n-1);           % Producers in compartment n
    F = y(2*n);             % Free-rider in compartment n
    
    % dydt = RHS
    if P + F ~= 0
        dydt(2*n-1,1) = G1(b*P/(P+F),P+F,beta)*P*(1-Y/K) - delta_P*P;
        dydt(2*n,1)   = G2(b*P/(P+F),P+F,beta)*F*(1-Y/K) - delta_F*F;
    else
        dydt(2*n-1,1) = 0;
        dydt(2*n,1)   = 0;
    end

end

end

function number_in_compartment = mnrnd(Ncomp,x,prob)

 F = cumsum(prob);
 
 number_in_compartment = zeros(Ncomp,1);
 
 for m = 1 : x
    roll = rand;
    ind = find(roll<F,1);
    number_in_compartment(ind) = number_in_compartment(ind) + 1;
 end

end