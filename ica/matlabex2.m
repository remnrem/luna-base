POINTS = 1000; % number of points to plot

% define the two random variables
% -------------------------------
for i=1:POINTS
	A(i) = round(rand*99)-50;              % A
	B(i) = round(rand*99)-50;              % B
end;
figure; plot(A,B, '.');                        % plot the variables
set(gca, 'xlim', [-80 80], 'ylim', [-80 80]);  % redefines limits of the graph

% mix linearly these two variables
% --------------------------------
M1 = 0.54*A - 0.84*B;                          % mixing 1
M2 = 0.42*A + 0.27*B;                          % mixing 2
figure; plot(M1,M2, '.');                      % plot the mixing
set(gca, 'ylim', get(gca, 'xlim'));            % redefines limits of the graph

% withen the data
% ---------------
x = [M1;M2];
c=cov(x')		 % covariance
sq=inv(sqrtm(c));        % inverse of square root
mx=mean(x');             % mean
xx=x-mx'*ones(1,POINTS); % subtract the mean
xx=2*sq*xx;              
cov(xx')                 % the covariance is now a diagonal matrix
figure; plot(xx(1,:), xx(2,:), '.');

% show projections
% ----------------
figure; 
axes('position', [0.2 0.2 0.8 0.8]); plot(xx(1,:), xx(2,:), '.'); hold on;
axes('position', [0   0.2 0.2 0.8]); hist(xx(1,:));  set(gca, 'view', [90 90]);
axes('position', [0.2 0   0.8 0.2]); hist(xx(2,:));

% show projections
% ----------------
figure; 
axes('position', [0.2 0.2 0.8 0.8]); plot(A,B, '.'); hold on;
axes('position', [0   0.2 0.2 0.8]); hist(B);
set(gca, 'xdir', 'reverse'); set(gca, 'view', [90 90]);
axes('position', [0.2 0   0.8 0.2]); hist(A);

