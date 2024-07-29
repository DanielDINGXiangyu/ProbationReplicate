function [] = solve_nocorr(gamma,psi,filename)

%% Slove Bond Policy Function: Sudden stop with trend shock 
% Xiangyu DING (dingxiangyu@connect.hku.hk)


%% 0. Housekeeping
delete(gcp('nocreate')) % clear Parpoll

if verLessThan('matlab','8.2')
	% Old version. Use matlabpool
	matlabpool open
else
	% Newer version, can use parpool
	poolpbject=parpool(14);
end


%% 1. Values of model and algorithm parameters
%% 1.1 Epstein-Zin Utility
% gamma=10; %coefficient of relative risk aversion 
% psi=1.5; %coefficient of intertemporal elasticity of substitution

%% 1.2 Stochastic Process and grid

% 1.2.1 transitory relative tradable nontradable shock
StoVolZ=true; % true stochastic volatility, false constant volatility
rho_z_level=0.7501; % autocorrelation of transitory tradeable endowment % Seoane and Yuirdagul (2019) JIE
sig_z_uncon=0.0532; % UNCONDITIONAL standard deviation of tradeable endowment SHOCK (u_t) (not endowment itself) % Seoane and Yuirdagul (2019) JIE
rho_z_vol=0.90^4; % autocorrelation of the tradeable volatility % Colacito et al. 2022 RFS
sig_z_vol=0.15*2; % UNCONDITIONAL standard deviation of tradeable volatility SHOCK (u_t) % Colacito et al. 2022 RFS

% 1.2.1 growth rate shock (short run risk or trend)
StoVolG=true; % true stochastic volatility, false constant volatility
rho_g_level=0.5499; %autocorrelation of persistent growth trend G % Seoane and Yuirdagul (2019) JIE
sig_g_uncon=0.0353; % UNCONDITIONAL standard deviation of persistent growth trend G SHOCK (not endowment itself) % Seoane and Yuirdagul (2019) JIE
rho_g_vol=0.90^4; % autocorrelation of the growth volatility % Colacito et al. 2022 RFS
sig_g_vol=0.15*2; % UNCONDITIONAL standard deviation of growth volatility SHOCK (u_t) % Colacito et al. 2022 RFS

% StoR=false; % add stochastic interest rate in future

% The process is as follows:
% Y^T_t = Gamma_t Z_t; (Z_t is relative tradable cycle shock)
% Y^N_t = Z_t;
% Gamma_t/Gamma_{t-1} = G_t  (this is trend shock)

% log(x_t) = lambda*log(x_{t-1}) + u_t       where x \in {Gamma, Z} 
% vol_t = (1-rho)*mue + rhoe*vol_{t-1} + epsilon_t
% u_t ~ N(0,exp(vol_t)); epsilon_t ~ N(0,sigmae^2)
%% 1.2 Other Parameter 
ita=(1/0.83)-1; %CES preferences exponent parameter (elasticity=1/(1+ita))
beta=0.91; %discount factor
omega=0.31; % CES preferences share parameter
kappa=0.32; %fraction of income pledgeable as collateral
r=0.04; % world interest rate

%% 1.3 Grid for exogeous state and endogenous state
% 1.3.1 exogenous state
NZlevel=5;
NZsigma=3;
NGlevel=5;
NGsigma=3;
% 1.3.2 endogenous state: bond holding
NB=200;
%% 1.4 Algorithm parameters

% 1.4.1 Iteration parameters
% Set uptd, updt_SPP small if not converge.
% Set VFI=true for VFI in each loop (check convergece)

VFI=false; % if VFI=true, then value function iteration inside policy iteration loop (slow but robust)
uptd=0.2; %0.1; %Weight on new policy function in update to next iteration
updt_SPP=0.005;  %Weight on new policy function: Social planner equilibrium

outfreq=20; %Display frequency (shows in screen each 20th iteration)
iter_tol=10000;  % Maximum number of fixed-point iterations
tol=1e-6; %1e-4 1e-7; % Convergence criterion (lower for session, higher for accuracy
tol_EEbind = 1.0e-16; % tolerance for binding constriant 


% 1.4.2 Ad-hoc borrowing constaint
bmin=-1.00; %-1.06; 
bmax=0.00; %-0.41;

%% 1.5 Define functional form to simplify code
% epstein-zin value function  fV(c_total,expec_V,growth)
fV = @(c_total,expec_V,growth) ((1-beta).*c_total.^(1-1/psi) + beta.*expec_V.^(1-1/psi).*growth.^(1-1/psi)).^(1/(1-1/psi));
% consumption bundle  fC(c_trade,c_nontrade)
fC = @(c_trade,c_nontrade) (omega.*c_trade.^(-ita) + (1-omega).*c_nontrade.^(-ita)).^(-1/ita);
% marginal utility of consumption fMup(c_trade,c_nontrade)
fMup = @(c_trade,c_nontrade) ((omega.*c_trade.^(-ita) + (1-omega).*c_nontrade.^(-ita)).^(-1/ita)).^(1+ita-1/psi).*omega.*c_trade.^(-ita-1);
% pricing fP(c_trade,c_nontrade)
fP = @(c_trade,c_nontrade)  (1-omega)./omega.*(c_trade./c_nontrade).^(1+ita);


%% 2 Discretization
%%%%%%%%%%  2.1 Farmer & Toda (2017, Quantitative Economics) %%%%%%%%%%%%%
% add path
addpath(genpath("discretization")) ;

%% 2.1 stochastic v.s. constant growth volatility
if StoVolG==true           % stochastic volaility growth
    [Prob_growth,gvGrids] = discreteSV(rho_g_level,rho_g_vol,sig_g_uncon,sig_g_vol,NGlevel,NGsigma);
    NG = NGlevel*NGsigma;
else                       % constant volaility growth
    [gvGrids,Prob_growth] = tauchenhussey(NGlevel,0,rho_g_level,sig_g_uncon,sig_g_uncon);
    gvGrids=gvGrids';
    NG = NGlevel;
end

%% 2.2 stochastic v.s. constant tradable level volatility
if StoVolZ==true
    [Prob_trade,zvGrids] = discreteSV(rho_z_level,rho_z_vol,sig_z_uncon,sig_z_vol,NZlevel,NZsigma);
    NZ = NZlevel*NZsigma;
else
    [zvGrids,Prob_trade] = tauchenhussey(NZlevel,0,rho_z_level,sig_z_uncon,sig_z_uncon);
    zvGrids=zvGrids';
    NZ = NZlevel;
end

%% 2.3 Generate mat
NSS = NZ*NG; % Total Number of Exogenous States (s)

% 2.3.1 Generate exogenous state mat: z z_vol g g_vol
% T1G1, T1G2,...| T2G1, T2G2,| % This form is easy to compare the growth grid
zvGrids = kron(zvGrids, ones(1,NG)); % First row is level second row is std
gvGrids = kron(ones(1,NZ), gvGrids); % First row is level second row is std

% Large Transition Prob
Prob = kron(Prob_trade,Prob_growth); % Transition matrix of the

% Exogenous state_mat
z=exp(repmat(zvGrids(1,:),NB,1)); % output state 
g=exp(repmat(gvGrids(1,:),NB,1)); % growth state

yt=z.*g;
yn=g;


% interest rate mat
R=repmat(r,NB,NSS); % Interest Rate at each (b,z)
Q=1./(1+R);         % Bond price at each state (b,z)


% 2.3.2 Generate endogenous state mat: bond
% Create bonds grid
B=zeros(NB,1);
B(1)=bmin;
%%%%%%%%%%%%%%% Symmetric  %%%%%%%%%%%%%%%%%%%%%
%%%%%%% (finer grid for low i to improve accuracy near lower bound)%%%%%
for i=2:NB
    B(i)=B(i-1)+(bmax-B(i-1))/((NB-i+1)^1.05); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b=repmat(B,1,NSS);  % Debt at each state (b,z)
%% 3. Decentralized Equilibrium
%%  3.1 Initial values 
%%%%%%%%%%%%%%%%%%%  Step 1 of Cell 3 in the handout
%%%%%%%%%%Initial conjectures of decision rules 
%%%%%%(they are updated before the end of each iteration). 
bp=b;    % Initial bonds dec rule: b'(b,z)=b (45degree line)
ct=b+yt-Q.*bp.*g; % Initial cT dec rule implied by initial bonds dec rule
price=fP(ct,yn);   %(1-omega)/omega*(ct./yn).^(1+ita); % Initial Nontradables price implied by cT
bpmax=-kappa*(price.*yn +yt)./Q./g; %-kappa*(price+yt);  % Debt Limit
ctbind=b+yt-bpmax.*Q; % Initial constrained cT values

EE=zeros(NB,NSS);       % Initial Unconstrained Euler Equation Error
mup=ones(NB,NSS);       % Preallocation for marginal utility of tradables
emu=zeros(NB,NSS);      % Preallocation for expected marginal utility
V=ones(NB,NSS);
EV=zeros(NB,NSS); 

%% 3.2 Fixed point Iteration Loop (see Cell 3, Step 2 in Handout)
iter=0;  % Initialize iteration counter
dis=100;  % Initialize value of convergence metric (just to enter loop)
options =optimset('Display','off'); %options for fzero nonlinear solver
tic
disp('DE Iter      Norm');

     %%%%%%%%Step 2 of Cell 3 in handout: 
while dis>tol && iter <iter_tol  %Exits loop if either convergence is attained
    % or total number of iterations exceed value of iter_tol parameter
    
%    oldp=price; %conjectured pN for current iteration
    oldct=ct; % conjectures cT for current iteration
    oldbp=bp; %conjecture bonds dec rule for current iteration
    oldV=V;
    
    %% 3.4 Update Value function 
    % Flow utilities for each regime
    % Initialize expected utility at current flow utility
    
    % Two methods 
    if VFI~=true
        %% Method1 (fast): Update value function each policy iteration
        V_sigma=V.^(1-gamma);
        %totalc=fC(ct,yn); %(omega.*ct.^(-ita) + (1-omega).*yn.^(-ita)).^(-1/ita); 
        parfor i=1:NB
            for j=1:NSS
                % Calculate expected utility by interpolating current Value
                % function.
                EV(i,j)=(interp1(B,V_sigma,bp(i,j),'linear','extrap')*Prob(j,:)')...
                        .^(1/(1-gamma));
            end
        end
        V=fV(fC(ct,yn),EV,g); %((1-beta).*totalc.^(1-1/psi) + beta.*EV.^(1-1/psi).*g.^(1-1/psi)).^(1/(1-1/psi));
             
        
    else
        %% Method2 (slow but robust): Update Value function using value function iteration
        % Initiallize matrices for value functions of the Decentralize Equilibrium

        % Flow utilities for each regime
        %totalc=fC(ct,yn); %(omega.*ct.^(-ita) + (1-omega).*yn.^(-ita)).^(-1/ita); 

        dis2=100;                         % Initial convergence criterion
        iter2=0;                         % Iteration counter
        % disp('Value Iter     Norm');    
        while dis2>tol && iter2 <iter_tol
            V_sigma = V.^(1-gamma);
            parfor i=1:NB
                for j=1:NSS
                    % Calculate expected utility by interpolating current Value
                    % function.
                    EV(i,j)=(interp1(B,V_sigma,bp(i,j),'linear','extrap')*Prob(j,:)')...
                        .^(1/(1-gamma));
                end
            end

            % New value function is flow utility plus NEW expected utility
            V_new=fV(fC(ct,yn),EV,g); % ((1-beta).*totalc.^(1-1/psi) + beta.*EV.^(1-1/psi).*g.^(1-1/psi)).^(1/(1-1/psi));

            % Check convergence for value functions
            dis2=max(max(abs((V_new-V)./V)));

            iter2=iter2+1;                    % Update counters
            V=V_new;                        % Uptade Value Functions
            % if mod(iter2, outfreq) == 0;
            %     fprintf('%d          %1.7f \n',iter2,d3); % Display every (outfreq) iterations
            % end

        end
        % fprintf('%d          %1.7f \n',iter2,d3); % Display last iteration
    end

    %% 3.5 Calculate the LHS of Euler Equation
    
    % totalc=fC(ct,yn); %(omega.*ct.^(-ita)+(1-omega).*yn.^(-ita)).^(-1/ita);   
    % mup(b,z) is Marginal utility of C^T given current consumption policy c(b,z)
    mup=fMup(ct,yn); %totalc.^(1+ita-1/psi).*omega.*(ct.^(-ita-1)); 
    
    % V_sigma is the term inside exptctation operator of the expectation term
    V_sigma=V.^(1-gamma);
    
    % Co is the term inside the expectation operator of the covariance term
    Co=V.^(1/psi-gamma).*mup;
    parfor i=1:NB % Can be parallelized using for instear of for
        for j=1:NSS                
            %EMU is expected marginal utility tomorrow in today's grid, given
            %current future debt policy bp(b,z). Found by interpolation
            emu(i,j)=beta*(1+R(i,j))...
                *(interp1(B,V_sigma,bp(i,j),'linear','extrap')*Prob(j,:)')...
                    ^((gamma-1/psi)/(1-gamma))... % Expectation term
                *(interp1(B,Co,bp(i,j),'linear','extrap')*Prob(j,:)')... % covariance term
                *g(i,j)^(-1/psi); % current grwoth term

%             CE = (interp1(B,V_sigma,bp(i,j),'linear','extrap')*Prob(j,:)')...
%                     ^((gamma-1/psi)/(1-gamma));
%             emu(i,j)=beta*(1+R(i,j))*...
%                 (interp1(B,Co./CE,bp(i,j),'linear','extrap')*Prob(j,:)'); 
        end
    end
    
    %% 3.6 Calculate RHS of Euler, and the Euler function gap
    parfor i=1:NB % % Can be parallelized using for instear of for
        for j=1:NSS  
      %%%%%%%%%%%% Step 3 of Cell 3 in the handout.
          % Calculate Euler Equation gap assuming the constraint binds
          % "today" (i.e. using cbind(b,z) for uT(t), this is eq. (4) in the handout. 
            EE(i,j)=fMup(ctbind(i,j),yn(i,j))-emu(i,j); 
            %(omega*ctbind(i,j)^(-ita)+(1-omega)*yn(i,j)^(-ita))^(1/(ita*psi)-1/ita-1)*omega*ctbind(i,j)^(-ita-1)-emu(i,j);
           
            
           %% 3.6.1 binding state
            % first term of EE: fMup(fC(ctbind(i,j),yn(i,j)),ctbind(i,j))
            % is the lowest possible marginal utility of tradable consumption, by consumption up to constraint
            % if still larger than emu (highest consumption cannot equate the euler), then constaint must binding
            
            cinit=ct(i,j); %initial condition for solver  
            if EE(i,j)>tol_EEbind %tol*1e-5	% If positive, constraint will be binding then:
                
                bp(i,j)=bpmax(i,j);  % bonds determined by binding credit constraint
                ct(i,j)=ctbind(i,j);   % consumption solved from budget constraint
          %% 3.6.2 Slack state
            else % Constraint not binding
                % Define function used to find consumption decision
                % that solves Euler Equation 
                ff = @(cc) fMup(cc,yn(i,j))-emu(i,j);  
                % Solve Euler Equation to get new consumption dec rule
                c0=cinit;
                [ct(i,j),EE(i,j)]=fzero(ff,c0,options); 
                % Solve bonds dec rule from budget constraint (and within
                % grid boundaries) 
                bp(i,j)=max((1/Q(i,j)/g(i,j))*(yt(i,j)+b(i,j)-ct(i,j)),bmin);
                bp(i,j)=min(bp(i,j),bmax); 
                % Here debt may not be consistent with consumption anymore
                % but it will be the optimal choice given the grid.       
            end
        end
    end
    price=fP(ct,yn); % (1-omega)/omega*(ct./yn).^(1+ita);  % set PN to clear the N market
    
    %% 3.7 Update Policy function
    %===============Consolidate new solutions==============================
    % To make sure c is feasible:
    ct=b+yt-max(Q.*bp.*g,-kappa*(yt+price.*yn));    % qb'>=-kappa*(yt+yn*pn)
    price=fP(ct,yn); %(1-omega)/omega*(ct./yn).^(1+ita); % Price consistent with c
    bp=(1./Q./g).*(b+yt-ct); % In principle this line should be redundant
    % End of Steps 3,4,5 of Cell 3 of the handout
    %========================================================================
 
    %%%Step 6 of Cell 3 of handout
    
    %=====================Updating rules for next iteration==============
    bp=uptd*bp+(1-uptd)*oldbp;
    ct=uptd*ct+(1-uptd)*oldct;
%    price=uptd*price+(1-uptd)*oldp;
    %========================================================================    
 
   %%%%%%%%%%%%%% Find new constrained values:%%%%%%%%%%%%%%%%%%%
    %%%%% This updates constrained dec rules for next iteration.
    bpmax=-kappa*(yt+price.*yn)./Q./g; %-(1./Q).*kappa.*(price+yt);
    bpmax(bpmax>bmax)=bmax;
    bpmax(bpmax<bmin)=bmin;
    ctbind=b+yt-Q.*bpmax.*g;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    iter=iter+1; % Update iteration counter
    
    %Calculate difference between new and old policies for convergence
    %check
    dis=max([max(max(abs(ct-oldct))),max(max(abs(bp-oldbp))),max(max(abs(V-oldV)))]);
    
    % Print results once every (outfreq) iterations
    if mod(iter, outfreq) == 0;
        fprintf('%d          %1.7f \n',iter,dis);
    end
end
fprintf('%d          %1.7f \n',iter,dis);    % Print last iteration result
toc

%% ---------------------------------------------
%% 4. Social Planner's Solution (SP) (Cell 4)
uptd=updt_SPP;    % Slower updating speed of dec rule for planners problem 
%% 4.1 Initial Conjecture 
%%%%%%%Initial conjectures for SP solution set at DE solutions or same as
%%%%%%%in DE case
bpSP=bp; %bp or b
ctSP=b+yt-Q.*bpSP.*g;  %c;
priceSP=fP(ctSP,yn);    %(1-omega)/omega*(ctSP./yn).^(1+ita); %price;

bpmaxSP=-kappa*(priceSP.*yn +yt)./Q./g;  %bpSP;
cbindSP=b+yt-Q.*bpmaxSP.*g;  %c;
EESP=EE;

mupconsSP=ones(NB,NSS); % Preallocation for SP's marginal utility of tradables from consumption
mupSP=ones(NB,NSS); % Preallocation for SP's marginal utility of tradables
emuSP=zeros(NB,NSS); % Preallocation for SP's exp. marginal utility

VSP=V;
EVSP=EV; 
%%%%%%%


%% 4.2 Policy function iteration
dis=100;      % Initialize difference between policies (just to enter loop)
iter=0;      % Initialize Iteration Counter
disp('SPP Iter     Norm');

while dis>tol && iter <iter_tol*2 %Exit if either convergence is attained.
    %or loop runs for more than number of times set in iter_tol parameter
    % times 2. Uses double because of slower updating
    
    % Set conjectures for current iteration
%    oldpSP=priceSP;
    oldcSP=ctSP;
    oldbpSP=bpSP;
    oldVSP=VSP;
	
    
    
    %% 4.3 Update Value function
     if VFI~=true
        
        %% Method1 (fast): Update value function each policy iteration
        VSP_sigma = VSP.^(1-gamma);
        totalcSP=fC(ctSP,yn); 
        parfor i=1:NB
            for j=1:NSS
                % Calculate expected utility by interpolating current Value
                % function.
                EVSP(i,j)=(interp1(B,VSP_sigma,bpSP(i,j),'linear','extrap')*Prob(j,:)')...
                        .^(1/(1-gamma));
            end
        end
        VSP=fV(totalcSP,EVSP,g);
        
    else
        %% Method2 (slow but robust): Update Value function using value function iteration
        % Initiallize matrices for value functions of the Decentralize Equilibrium

        % Flow utilities for each regime
        totalcSP=fC(ctSP,yn); 

        dis2=100;                         % Initial convergence criterion
        iter2=0;                         % Iteration counter
        % disp('Value Iter     Norm');    
        while dis2>tol && iter2 <iter_tol
            VSP_sigma = VSP.^(1-gamma);
            parfor i=1:NB
                for j=1:NSS
                    % Calculate expected utility by interpolating current Value
                    % function.
                    EVSP(i,j)=(interp1(B,VSP_sigma,bpSP(i,j),'linear','extrap')*Prob(j,:)')...
                        .^(1/(1-gamma));
                end
            end

            % New value function is flow utility plus NEW expected utility
            VSP_new=fV(totalcSP,EVSP,g);

            % Check convergence for value functions
            dis2=max(max(abs((VSP_new-VSP)./VSP)));

            iter2=iter2+1;                    % Update counters
            VSP=VSP_new;                        % Uptade Value Functions
            % if mod(iter2, outfreq) == 0;
            %     fprintf('%d          %1.7f \n',iter2,d3); % Display every (outfreq) iterations
            % end

        end
        % fprintf('%d          %1.7f \n',iter2,d3); % Display last iteration
    end
	
    % Find marginal utility of tradables consumption

    mupconsSP=fMup(ctSP,yn);
	
    % Calculate externality term, Equation (8) in the handout:
    PsiSP=kappa*(1-omega)/omega*(ita+1)*(ctSP./yn).^ita;
    PsiSPbind=kappa*(1-omega)/omega*(ita+1)*(cbindSP./yn).^ita;
    % Calculate social marginal utility, equation (7) in the handout
    mupSP=EESP.*PsiSP+mupconsSP;
    
    % VSP_sigma is the term inside exptctation operator of the expectation term
    VSP_sigma=VSP.^(1-gamma);
    
    % CoSP is the term inside the expectation operator of the covariance term
    CoSP=VSP.^(1/psi-gamma).*mupSP;
    
    %%%% Find expected marginal utility tomorrow in today's grid, given
    %%%% current future debt policy bp(b,z). By interpolation
    parfor i=1:NB % % Can be parallelized using for instear of for
        for j=1:NSS
            emuSP(i,j)=beta*(1+R(i,j))...
                *(interp1(B,VSP_sigma,bpSP(i,j),'linear','extrap')*Prob(j,:)')...
                    ^((gamma-1/psi)/(1-gamma))... % Expectation term
                *(interp1(B,CoSP,bpSP(i,j),'linear','extrap')*Prob(j,:)')...% covariance term
                *g(i,j)^(-1/psi); % current grwoth term
			%EMU is expected marginal utility tomorrow as function of (b,z)
        end
    end
    
    parfor i=1:NB % % Can be parallelized using for instear of for
        for j=1:NSS
            
            %%%%%%% Step 2: Calculate Euler equation error at max debt.
            EESP(i,j)=(fMup(cbindSP(i,j),yn(i,j))-emuSP(i,j))/(1-PsiSPbind(i,j));
            
            
            if EESP(i,j)>tol_EEbind % Constraint binds: go to step 3
                
                %%%% Step 3
                bpSP(i,j)=bpmaxSP(i,j); % debt will be as big as possible
                ctSP(i,j)=cbindSP(i,j);  % consumption is solved from budget constraint
               
            else % Debt restriction does not bind: go to step 4 in handout
                
                %%%%% Step 4 ,section 3 of handout:
                
                % Define function used to find consumption decision
                % that solves planner's Euler Equation 
                ff = @(cc) fMup(cc,yn(i,j))-emuSP(i,j); 
                % Solve consumption from Euler Equation
                cinitSP=ctSP(i,j); %initial condition for solver
                c0=cinitSP;
                [ctSP(i,j),EESP(i,j)]=fzero(ff,c0,options); 
                % Find planner's bond decision and make sure is inside the grid
                bpSP(i,j)=max((1/Q(i,j)/g(i,j))*(yt(i,j)+b(i,j)-ctSP(i,j)),bmin);
                bpSP(i,j)=min(bpSP(i,j),bmax);
            end
        end
    end
    priceSP= fP(ctSP,yn); 
    
    %===============Check collateral constraint==============================
    % Again debt was forced into the grid but may not be consistent with
    % consumption.
    ctSP=b+yt-max(Q.*bpSP.*g,-kappa*(yt+priceSP.*yn));
    priceSP= fP(ctSP,yn); 
    bpSP=(b+yt-ctSP)./Q./g;
    %========================================================================
    
    %=====================Updating rule.Must be slow, important!==============

%if d2<1e-2 % Smaller Updates as dec rules are closer to converging 
 %   uptd=0.001;
%else
 %   uptd=0.005;
%end  
   
    %=====================Updating rule==============
    bpSP=uptd*bpSP+(1-uptd)*oldbpSP;
    ctSP=uptd*ctSP+(1-uptd)*oldcSP;
%    priceSP=uptd*priceSP+(1-uptd)*oldpSP;
    %========================================================================
 
       %%%%%%%%% Find constrained values
     bpmaxSP=-(1./Q./g).*kappa.*(priceSP.*yn+yt);
     bpmaxSP(bpmaxSP>bmax)=bmax;
     bpmaxSP(bpmaxSP<bmin)=bmin;
     cbindSP=b+yt-Q.*bpmaxSP.*g;

    %%%%%%%%%%%%%  
        
    % Determine distance between policies
    dis=max([max(max(abs(ctSP-oldcSP))),max(max(abs(bpSP-oldbpSP)))]);
    
    iter=iter+1; %update iteration counter
    iter_sp = iter;
    d2_sp   = dis;
    if mod(iter, outfreq*5) == 0; %display every 5*outfreq results
        %fprintf('%d          %1.7f \n',iter,d2);
        fprintf('%d          %1.7f \n',iter,d2_sp);
        %fprintf('%d          %10.2e \n',iter_sp,d2_sp);
    end   
end


%% 5. Alternative way of detrend: x_curr=x/current trend
V_curr = V./g;
VSP_curr = VSP./g;
b_curr = b./g;

%% 6. Social Welfare Gain
Wf=(VSP./V)-1;

%% 7. Optimal Tax

% Notice that mupSPtax different from mupSp: It does not add EESP*psi
mupSPtax=fMup(ctSP, yn);  %mup(b,y)

% VSP_sigma is the term inside exptctation operator of the expectation term
VSP_sigma_tax=VSP.^(1-gamma);
    
% CoSP is the term inside the expectation operator of the covariance term
CoSP_tax=VSP.^(1/psi-gamma).*mupSPtax;

%%%% Find expected marginal utility tomorrow in today's grid, given
%%%% current future debt policy bp(b,z). By interpolation
parfor i=1:NB % % Can be parallelized using for instear of for
    for j=1:NSS
        emuSP_tax(i,j)=beta*(1+R(i,j))...
            *(interp1(B,VSP_sigma_tax,bpSP(i,j),'linear','extrap')*Prob(j,:)')...
                ^((gamma-1/psi)/(1-gamma))... % Expectation term
            *(interp1(B,CoSP_tax,bpSP(i,j),'linear','extrap')*Prob(j,:)')... % covariance term
            *g(i,j)^(-1/psi);
        %EMU is expected marginal utility tomorrow as function of (b,z)
    end
end

% Optimal tax: find the gross tax rate, equation (11) in the handout
Tau=mupSP./emuSP_tax-1;

Tau(EESP>0.001)=0;
Tau(Tau<0)=0;

%% 7.2 Another way of solving the optimal tax

CoSP_gap = VSP.^(1/psi-gamma).*EESP.*PsiSP; % numerator 
CoSP_tax = VSP.^(1/psi-gamma).*mupSPtax;    % denomenitor
parfor i=1:NB % % Can be parallelized using for instear of for
    for j=1:NSS
        Tau_num(i,j) = (interp1(B,CoSP_gap,bpSP(i,j),'linear','extrap')*Prob(j,:)');
        Tau_den(i,j)  = (interp1(B,CoSP_tax,bpSP(i,j),'linear','extrap')*Prob(j,:)'); 
    end
end
Tau_2 = Tau_num./Tau_den;

%figure;
%plot(B,Tau); legend('1','2','3','4','5','6','7','8','9')
%xlabel('debt')

 %title('tax '); ylabel('%')
% Recall this is tax only for non-binding region for tax
% When unconstrained future debt increases with current debt
% Tax is set to 0 when constraint binds

AN=bpSP(2:end,:)-bpSP(1:end-1,:);   % change in future debt for adjacent current debt

for i=1:NSS
    % Find the first positive change, the first unconstrained debt choice
    IN(i)=find(AN(:,i)>0,1,'first'); 
end


for j=1:NSS
    % For all constrained, set tax to 1
    if IN(j)>1
        Tau(1:IN(j)-1,j)=0;
        Tau_2(1:IN(j)-1,j)=0;
    end
end

% just in case, all taxes should be 0 or more.
Tau(Tau<0)=0;
Tau_2(Tau_2<0)=0;


%% 8. Save Result
save(append('Solution/LowBound/nocorr_',filename))  
toc

end