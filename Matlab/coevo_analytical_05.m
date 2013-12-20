% this is an attempt to make the model analytical. It is assumed that the
% performance of individual learners is constant. Path dependence is
% ignored. Social learners have 3 individuals as sample. Choice proportions
% are always in equilibrium.

clear all

% PARAMETERS
tmax=50;                    % periods per generation
gen=200000;                 % generations
selCoef=1/10;               % maximum selection coefficient
nStrats=7;                  % total possible number of strategies
vertTransm=1;               % is there vertical tranmission from parents to offspring?
regime=1;                   % how the environment changes. 0->no regression to the mean, 1->medium, 32->high
incr=2/100;                 % increment at which the environment becomes better or worse
pincr=1;                    % probability that environmental quality changes at all after each period
tic

% STRATEGIES
% 1 ind: individual
% 2 con: conformists
% 3 pbsl1: payoff-biased social learners equal weight
% 4 pbsl3: payoff-biased social learners, more weight on gains
% 5 mce: McElreath payoff-biased social learners
% 6 OIL: opportunistic individual learners
% 7 OC: opportunistic conformist

xInit=zeros(nStrats,1);

% xInit(1)=.5;
% xInit(2)=.5;
% xInit(3)=.0;
% xInit(4)=.0;
% xInit(5)=.0;
% xInit(6)=.0;
% xInit(7)=.0;

% xInit(1)=1/2;
% xInit(2)=1/2;
% xInit(3)=.0;
% xInit(4)=.0;
% xInit(5)=.0;
% xInit(6)=.0;
% xInit(7)=.0;

xInit=1/nStrats;

% TESTED PARAMETER RANGE
incrVec=[.01:.005:.1];;
pincrVec=[.1:.05:1];
% incrVec=.02;
% pincrVec=1;
% tmaxVec=[10:10:150];
% tmaxVec=50;

% define the FOCAL PARAMETER VECTORS
u2=incrVec;
u1=pincrVec;
% u1=nRoundVec;
ylab='p_{incr}';
xlab='k_{incr}';

nRepit=length(u1)*length(u2);

m=zeros(nStrats,length(u1)*length(u2));
m=reshape(m,nStrats,length(u1),length(u2));
mv=m;   %variance

wb=waitbar(0,'progress');

for j=1:length(u1)
    
    pincr=pincrVec(j);
    
    for i=1:length(u2)
        
        incr=incrVec(i);
        
        pA0=1/2;pB0=1/2+incr/2;     % initial values of pA and pB, so that pA=/=pB forall t
            
        [x c w maxGen pA pB]=coevo_ana_05(...
        xInit,...                   % initial conditions
        tmax,...                    % periods per generation
        gen,...                     % generations
        selCoef,...                 % maximum selection coefficient
        nStrats,...                 % total possible number of strategies
        vertTransm,...              % whether there is vertical transmission
        regime,...                  % how the environment changes. 0->no regression to the mean, 1->medium, 32->high
        pA0,...
        pB0,...                     % initial values of pA and pB
        incr,...                    % increment at which the environment becomes better or worse
        pincr);                     % probability that environmental quality changes at all after each period

        xLast=mean(x(:,maxGen-9999:maxGen),2);
        xVar=var(x(:,maxGen-9999:maxGen)')';
        m(:,j,i)=xLast;
        mv(:,j,i)=xVar;
        
        if length(u1)>1
            waitbar(j/length(u1),wb);
        else
            waitbar(i/length(u2),wb);
        end
        
    end
end

% FIGURES


if nRepit>=6    
    
    for s=1:nStrats
        % the focal strategy
        mc=squeeze(m(s,:,:));
        figure
        % colormap from white->low freq to black->high freq
        colormap([1-[0:.01:1]' 1-[0:.01:1]' 1-[0:.01:1]']);
        % add 0 and 1 to data so that full scale [0;1] is used
        mc2=[mc;zeros(1,length(u2))];
        mc2(length(u1)+1,length(u2))=1;
        lastEntry=2*u1(length(u1))-u1(length(u1)-1);
        % plot "heat map"
        surf(u2,[u1 lastEntry],zeros(length(u1)+1,length(u2)),mc2);
        axis([min(u2) max(u2) min(u1) max(u1) 0 1])
        set(gca,'fontsize',14)
        xlabel(xlab,'fontsize',14)
        ylabel(ylab,'fontsize',14)
        % view from above
        view([0 90])
        axis square
        colorbar
        title(s)
    end
    
elseif nRepit==1
    figure
    hold on
    plot([0:maxGen],x')
    axis([0 maxGen 0 1])
    set(gca,'fontsize',14)
    
%     figure
%     scatter(x(4,2:gen),(w(4,2:gen)-1)/selCoef,'k.')
    
end

toc

close(wb)