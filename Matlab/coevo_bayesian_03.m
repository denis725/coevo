% this program tests bayesian learning in the coevo context
% in this version, at the beginning of each generation
% pA and pB are drawn randomly according to their probability distribution
% but remain constant thereafter during the whole generation.

clear all

% 'random seed'
% rand('seed',17413)

% PARAMETERS
tmax=500;
iterations=1;
nindi=1000;

% these are the parameters used to estimate the prob. distrib.
incr=0.02;
pincr=1;
regime=1;
dp=0;
pA0=0.5;
pB0=0.5;


% load prior probability distribution of pA and pB:
load('probdistp');
pAbIni=probdistp;
pBbIni=probdistp;
cprob=cumsum(probdistp);

% the single probabilities:
probs=[0:incr:1];


% INITIALIZE

% prior probabilities of success and failure:
pAs=sum(probs.*pAbIni');
pAf=sum((1-probs).*pAbIni');
pBs=sum(probs.*pBbIni');
pBf=sum((1-probs).*pBbIni');

% choices
nChoice=zeros(nindi,tmax*iterations);
% 1st choice random
nChoice(:,1)=(rand(nindi,1)>.5);
% successes
nSuccess=zeros(nindi,tmax*iterations);
% best choice made?
nBest=zeros(nindi,tmax*iterations);
% expectations: mean and median
expAmean=zeros(1,tmax*iterations);
expBmean=zeros(1,tmax*iterations);
expAmed=zeros(1,tmax*iterations);
expBmed=zeros(1,tmax*iterations);

% matrix containing performances
perfmat=zeros(nindi,tmax*iterations);

% total environment:
superpA=zeros(1,tmax*iterations);
superpB=zeros(1,tmax*iterations);

wb=waitbar(0,'progress');

for i=1:iterations
    
    % establish prior at beginning of each generation
    pAb=zeros(nindi,length(pAbIni),tmax);
    pBb=zeros(nindi,length(pBbIni),tmax);
    pAb(:,:,1)=repmat(pAbIni,1,nindi)';
    pBb(:,:,1)=repmat(pBbIni,1,nindi)';
    
    % generate environment
    randA=rand;randB=rand;
    rindA=find(cprob>randA);
    rindB=find(cprob>randB);
    pA=probs(rindA(1))*ones(1,tmax);
    pB=probs(rindB(1))*ones(1,tmax);
    
    %random matrix used for determining success
    randsucmatA=rand(nindi,tmax);
    randsucmatB=rand(nindi,tmax);
    
    for t=1:tmax
        
        % determine success
        nSuccess(:,t+(i-1)*tmax)=nChoice(:,t+(i-1)*tmax).*(randsucmatA(:,t)<pA(t))+...
            mod(nChoice(:,t+(i-1)*tmax)+1,2).*(randsucmatB(:,t)<pB(t));
        
        % UPDATE THE PRIORS
        % determine indices of the different the outcomes
        indAs=(nChoice(:,t+(i-1)*tmax)==1)&(nSuccess(:,t+(i-1)*tmax)==1);
        indAf=(nChoice(:,t+(i-1)*tmax)==1)&(nSuccess(:,t+(i-1)*tmax)==0);
        indBs=(nChoice(:,t+(i-1)*tmax)==0)&(nSuccess(:,t+(i-1)*tmax)==1);
        indBf=(nChoice(:,t+(i-1)*tmax)==0)&(nSuccess(:,t+(i-1)*tmax)==0);
        % new priors
        pAb(indAs,:,t+1)=pAb(indAs,:,t).*repmat(probs,sum(indAs),1)/pAs;
        pAb(indAf,:,t+1)=pAb(indAf,:,t).*(1-repmat(probs,sum(indAf),1))/pAf;
        pBb(indBs,:,t+1)=pBb(indBs,:,t).*repmat(probs,sum(indBs),1)/pBs;
        pBb(indBf,:,t+1)=pBb(indBf,:,t).*(1-repmat(probs,sum(indBf),1))/pBf;
        % priors of the non chosen option remain the same
        pBb(indAs,:,t+1)=pBb(indAs,:,t);
        pBb(indAf,:,t+1)=pBb(indAf,:,t);
        pAb(indBs,:,t+1)=pAb(indBs,:,t);
        pAb(indBf,:,t+1)=pAb(indBf,:,t);
        % normalize priors (should not be necessary, but matlab is too
        % imprecise)
        pAbsum=repmat(sum(pAb(:,:,t)'),length(probs),1);
        pBbsum=repmat(sum(pBb(:,:,t)'),length(probs),1);
        pAb(:,:,t)=pAb(:,:,t)./pAbsum';
        pBb(:,:,t)=pBb(:,:,t)./pBbsum';
        
        % decision: choose option with higher prior:
        nChoice(:,t+(i-1)*tmax+1)=(1+sign(sum(pAb(:,:,t)'.*repmat(probs,nindi,1)')'-sum(pBb(:,:,t)'.*repmat(probs,nindi,1)')'))/2;
        % if draw->random choice
        indDraw=(nChoice(:,t+(i-1)*tmax+1)==0.5);
        nChoice(indDraw,t+(i-1)*tmax+1)=(rand(sum(indDraw),1))>.5;
        
        % PERFORMANCE
        % determine whether best choice:
        dpApB=sign(pA(t)-pB(t));
        nBest(:,t+(i-1)*tmax)=((2*nChoice(:,t+(i-1)*tmax)-1)*dpApB);
        % performance of the strategy:
        perfmat(:,t+(i-1)*tmax)=(nBest(:,t+(i-1)*tmax)+1)/2;
        
        % DETERMINE WHAT IS BELIEVED TO BE THE TRUE VALUE
        expAmean(t+(i-1)*tmax)=mean(sum(pAb(:,:,t)'.*repmat(probs,nindi,1)')');
        expBmean(t+(i-1)*tmax)=mean(sum(pBb(:,:,t)'.*repmat(probs,nindi,1)')');
        expAmed(t+(i-1)*tmax)=median(sum(pAb(:,:,t)'.*repmat(probs,nindi,1)')');
        expBmed(t+(i-1)*tmax)=median(sum(pBb(:,:,t)'.*repmat(probs,nindi,1)')');
        
    end
    
    % remember last state
    pA0=pA(end);
    pB0=pB(end);
    
    % total environment:
    superpA(((i-1)*tmax+1):i*tmax)=pA;
    superpB(((i-1)*tmax+1):i*tmax)=pB;
    
    if i<iterations
        % inherit first choice from parent
        nChoice(:,1+tmax*i)=nChoice(:,tmax*i);
    end
    
    waitbar((i+1)/iterations,wb)
    
end

close(wb)

%plot part of the behavior
if iterations*tmax>=250
    tt=250;
    fig1=sum(nChoice)/nindi;
    fig1=fig1(1:end-1);
    figure
    plot([1:tt],fig1(1:tt),'k.')
    hold on
    plot([1:tt],superpA(1:tt)-superpB(1:tt)+.5,'k')
    plot([1 tt],[.5 .5],'k')
    axis([1 tt 0 1])
    set(gca,'fontsize',14)
    xlabel('period','fontsize',14)
    ylabel('proportion of A choices','fontsize',14)
end

% calculate performance
meanperf=mean(perfmat,1);
pident=find(superpA==superpB);
% remove draws
meanperf(pident)=[];
mean_performance=mean(meanperf)
std_error_mean=std(meanperf)*1.96/sqrt(iterations)


