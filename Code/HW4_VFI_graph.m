% HW4_VFI Cont. 

%Question 6 Plot Policy Function and Lorenz Curve
% Policy Function
figure(1);
plot(a,pol_fn);
title('Policy Function For each prductivity state');

% Gini coefficient and lorenz curve

popu = reshape(mu', [num_z * num_a, 1]);
wealth = reshape(repmat(a, num_z,1)',[num_z * num_a, 1]);

nonnegwealth = wealth;
nonnegwealth(nonnegwealth<0) = 0;

figure(2)
gini_wealth = gini(popu, nonnegwealth, true);
title({'Lorenz Curve for Wealth','Gini coefficient = ',num2str(gini_wealth)});

% wealth distribution
figure(3)
MU = sum(mu);
bar(a,MU)
title('Wealth Distribution')


% Track Euler Equation Error

Y = zgrid * w;
Y = repmat(Y',[1,num_a]);
A = repmat(a,[num_z,1]);
c=Y+ r * A - pol_fn;
cf=c(:,pol_indx');
cf1=reshape(cf,[num_z num_a num_z]);
i=1;
while i < num_z+1
c1(i,:)=PI(i,:)*cf1(:,:,i);
i=i+1;
end
Eulererror=sum(sum(abs(c.^(-sigma)-beta * c1.^(-sigma)* r ).*mu));

