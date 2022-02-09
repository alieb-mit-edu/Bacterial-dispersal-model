% competition simulations between two populations that are identical except
% for their detachment rates. 

%dt is detachment rate [hr-1]
%pop_1 is first population
%pop_2 is second population 
% *[] shows units

% the range of detachment rates 
dt = 0.01:0.05:1;
% the particle mortality rate (Mp) [hr-1]
Mp = 0.020;
 
%The mortality of free living cells [hr-1]
MF = 0.020;
 
% growth parameters: mumax is the maximum growth rate and B_threshold is
% the threshold number of cells on the particle that results in zero growth
% rate.
B_threshold = 5000000;
mumax = 0.45; %[hr-1]


% initial number of free living cells for each populations: [pop1 pop2]
BF_initial = [100000000 100000000];
 
% Number of particles 
particle_number = 1200;
 
%total number of time steps to run the simulations
T_total=100;
%time intervals [hr] 
T_interval = 0.5;
 
% particle and cell sizes
rc = 0.5*10^-6;
rp = 50*10^-6;

%cell-particle distance
rcp = 200*10^-6; %[m]
%%attachment probability [hr-1] calculated from the hitting probability 

%diffusion coefficient (cell and particle)
kB = 1.38*10^-23; % Boltzmann’s constant [J/K-1]
T = 293; % [K]
viscosity = 0.001; % [Pa s]
Dc = ((kB*T)/(6*pi*viscosity*rc))*3600;% [m2/hr]
Dp = ((kB*T)/(6*pi*viscosity*rp))*3600;% [m2/hr]

Pe=((rc+rp)./rcp).*erfc((rcp-(rc+rp))./sqrt(4.*(Dc+Dp).*T_interval));

% adjusting growth/mortality and attachment rates to given time intervals
mumax= T_interval.*mumax;
Pe=T_interval.*Pe;
Mp= T_interval.*Mp;
MF= T_interval.*MF;
 
for pop_2=1:length(dt)
    
    for pop_1=1:length(dt)
    
 % start the simulation from time step 1 
 time_step=1;
 
 % BF: number of free living cells at time step 1
 BF=BF_initial;
 
 % number of cells on particles are assumed to be zero at time zero. 
 BP=[];
 BP(1:particle_number,1:2)=0;
 
 % detachment rate of each population is expanded over individual particles
 % for easier matrix calculations
 dt_p=[dt(pop_1) dt(pop_2)].*ones(particle_number,2);
      
% adjusting the detachment rate to given time interval.
dt_p = T_interval.*dt_p;

 for time_step=1:1000

% mu is the actual growth rate calculated from negative liner relationship
% with the total number of cells on the particle
mu=mumax.*(1-sum(BP')./B_threshold); 
mu=[mu; mu]';

% if the number of cells on the particle exceeds a threshold value the
% growth rate becomes negative and we set the negative growth rate to zero
mu(mu<0)=0;

% Calculating the number of cells attached and detached to/from particles 

if Pe.*BF(1)<BF(1)/particle_number && Pe.*BF(2)<BF(2)/particle_number
    
 BP=BP+mu.*BP+Pe.*BF.*ones(particle_number,2)-dt_p.*BP;
 
 BF=BF-sum(Pe.*BF.*ones(particle_number,2))+sum(dt_p.*BP)-MF.*BF;
 
elseif Pe.*BF(1)>=BF(1)/particle_number && Pe.*BF(2)<BF(2)/particle_number
   
 
BP(:,1)=BP(:,1)+mu(:,1).*BP(:,1)+BF(1)/particle_number-dt_p(:,1).*BP(:,1);
 
 BF(1)=0;
 
  BP(:,2)=BP(:,2)+mu(:,2).*BP(:,2)+Pe.*BF(2).*ones(particle_number,1)-dt_p(:,2).*BP(:,2);
 
 BF(2)=BF(2)-sum(Pe.*BF(2).*ones(particle_number,1))+sum(dt_p(:,2).*BP(:,2))-MF.*BF(2);

 
 elseif Pe.*BF(1)<BF(1)/particle_number && Pe.*BF(2)>=BF(2)/particle_number
   
 
BP(:,2)=BP(:,2)+mu(:,2).*BP(:,2)+BF(2)/particle_number-dt_p(:,2).*BP(:,2);
 
 BF(2)=0;
 
  BP(:,1)=BP(:,1)+mu(:,1).*BP(:,1)+Pe.*BF(1).*ones(particle_number,1)-dt_p(:,1).*BP(:,1);
 
 BF(1)=BF(1)-sum(Pe.*BF(1).*ones(particle_number,1))+sum(dt_p(:,1).*BP(:,1))-MF.*BF(1);
    
 
   elseif Pe.*BF(1)>=BF(1)/particle_number && Pe.*BF(2)>=BF(2)/particle_number
   
 
     BP(:,2)=BP(:,2)+mu(:,2).*BP(:,2)+BF(2)/particle_number-dt_p(:,2).*BP(:,2);
 
 BF(2)=0;
 
   BP(:,1)=BP(:,1)+mu(:,1).*BP(:,1)+BF(1)/particle_number-dt_p(:,1).*BP(:,1);
 
 BF(1)=0;
 
end

% calculate particle mortality (particle_m). we uniform randomly choose from particles
% for predation. 
 particle_m = randsample(1:particle_number,round(particle_number.*Mp),true,sum(BP')./sum(sum(BP')));
 
% saving the total number of particle associated cells over time
BP_over_time(:,time_step)=[sum(BP(:,1)) sum(BP(:,2))];

% killing the particle associated cells of particles that are selected for predation. 
BP(particle_m,:)=0;

% plot population dynamics
 % figure(1)
 % semilogy(BP_over_time')
 % pause(0.0001)
 end

% relative abundance of population 1
Re_ab_pop_1(pop_1,pop_2)=BP_over_time(1,time_step)./sum(BP_over_time(:,time_step));

BP_over_time=[];
    end
   
end
imagesc(Re_ab_pop_1)
save Re_ab_pop_1 Re_ab_pop_1
 





