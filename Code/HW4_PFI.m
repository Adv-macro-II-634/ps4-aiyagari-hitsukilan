%%%%%HW3 Huggett's model %%%%%
clear;
clc;

% PARAMETERS
beta = .99; % discount factor 
sigma = 2; % coefficient of risk aversion
alpha = 1/3; % Cobb Douglas PF
l = 1; % unit of labor every work can supply
delta = .025; 
num_pfi = 30; % update policy fxn every k iteration

% Discretize exoge. st. varibale: labor productivity
num_z = 5;
rho = 0.5; % AR-1 
sigma_eps = 0.2; % AR-1 std. devs of residual
[zgrid, PI] = rouwenhorst(rho,sigma_eps,num_z);
zgrid = exp(zgrid);

PI_inv = PI^1000; %invariant distribution

N_s = zgrid * PI_inv(1,:)';

% Discretize endo.st.variable: asset
a_lo = 0; % no borrowing
a_hi = 100;% upper bound of grid points
num_a = 500;

a = linspace(a_lo, a_hi, num_a); % asset (row) vector

% INITIAL GUESS FOR aggregate K
K_min = 10;
K_max = 100;


K_tol= 1 ;

while (abs(K_tol) > 0.01)
    
    K_guess = (K_min + K_max) / 2;
    r = alpha * K_guess^(alpha-1) * N_s^(1-alpha) + (1-delta); %interest rate
    w = (1-alpha) * K_guess^alpha * N_s^(-alpha); %wage
    
    % CURRENT RETURN (UTILITY) FUNCTION
    
    cons = bsxfun(@minus, r * a', a);
    
    cons = bsxfun(@plus, cons, permute(zgrid, [1 3 2]) * w);
    
    ret = (cons .^ (1-sigma)) ./ (1 - sigma); % current period utility
    ret(cons<0) = -Inf;
    
    % INITIAL VALUE FUNCTION GUESS
    
    v_guess = zeros(num_z, num_a);
    
    % PFI
    v_tol = 1;
    
    while v_tol >.000001
        
        % CONSTRUCT RETURN + EXPECTED CONTINUATION VALUE
        v_mat = ret + beta * repmat(permute((PI * v_guess),[3 2 1]), [num_a 1 1]);

        % CHOOSE HIGHEST VALUE (ASSOCIATED WITH a' CHOICE)
       [vfn, pol_indx] = max(v_mat, [], 2); %max for each row
       
       v_tol = max(abs(permute(vfn, [3 1 2]) - v_guess));
       v_tol = max(v_tol(:));
       
       v_guess = permute(vfn, [3 1 2]);
       
       Q = makeQmatrix(pol_indx,PI);
       
       pol_fn = a(pol_indx);
       u_mat = bsxfun(@minus, r * a, pol_fn);
       u_mat = bsxfun(@plus, u_mat, zgrid' * w);
       u_mat = (u_mat.^(1-sigma)./(1-sigma));
       
       v_vec = v_guess(:);
       u_vec = u_mat(:);
       
       
       for j = 1:num_pfi
           up_v_vec = u_vec + beta * Q * v_vec;
           v_vec = up_v_vec;
       end
       
       v_guess = reshape(v_vec,num_z,num_a);
    
    end
    
    % KEEP DECSISION RULE
    pol_indx = permute(pol_indx, [3 1 2]);
   
    
    % SET UP INITITAL DISTRIBUTION
    mu = zeros(num_z,num_a);
    mu(:) = 1/(num_z * num_a);
    
    dis=1;
    
  while dis>0.0000001 
      % ITERATE OVER DISTRIBUTIONS
      MuNew = zeros(size(mu));
     [z_ind, a_ind, mass] = find(mu); % find non-zero indices
    
    for ii = 1:length(z_ind)
        apr_ind = pol_indx(z_ind(ii),a_ind(ii)); % which a prime does the policy fn prescribe?
        
        MuNew(:, apr_ind) = MuNew(:, apr_ind) + ... % which mass of households goes to which exogenous state?
            (P(z_ind(ii), :)*mass(ii))';
        
    end
    dis = max(max(abs(mu-MuNew)));
    mu = MuNew;
  end
  
  
   %Market clears
   aggK = sum(sum(mu.* pol_fn));
   K_tol = aggK-K_guess;
   if K_tol > 0
       K_min = K_guess;
   else
       K_max = K_guess;
   end
end






