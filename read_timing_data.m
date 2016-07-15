nmax = 9;

timing_coop_par.central = zeros(nmax,1);
timing_coop_par.init = zeros(nmax,2);
timing_coop_par.solve = zeros(nmax,2);
timing_coop_par.total = zeros(nmax,1);

timing_ncoop_par = timing_coop_par;

timing_coop_ser = timing_coop_par;
timing_ncoop_ser = timing_ncoop_par;

fname_par = 'parallel/timing_';
fname_ser = 'serial/timing_';

%% centralized
f = fopen([fname_par,'centralized.dat'],'r');
A = mean(fscanf(f,'%f',[3,Inf]),2);
fclose(f);
timing_centralized_par.central = A(1)*ones(nmax,1);
timing_centralized_par.init = A(2)*ones(nmax,1);
timing_centralized_par.solve = A(3)*ones(nmax,1);
timing_centralized_par.total = timing_centralized_par.central + ...
    timing_centralized_par.init + timing_centralized_par.solve;

f = fopen([fname_ser,'centralized.dat'],'r');
A = mean(fscanf(f,'%f',[3,Inf]),2);
fclose(f);
timing_centralized_ser.central = A(1)*ones(nmax,1);
timing_centralized_ser.init = A(2)*ones(nmax,1);
timing_centralized_ser.solve = A(3)*ones(nmax,1);
timing_centralized_ser.total = timing_centralized_ser.central + ...
    timing_centralized_ser.init + timing_centralized_ser.solve;

%% cooperative/non-cooperative

for i=1:nmax
    f = fopen([fname_par,'coop',num2str(i),'.dat'],'r');
    A = mean(fscanf(f,'%f',[5,Inf]),2);
    fclose(f);
    timing_coop_par.central = timing_coop_par.central + A(1);
    timing_coop_par.init = timing_coop_par.init + ones(nmax,1)*[A(2), A(4)];
    timing_coop_par.solve(i,:) = [A(3), A(5)];
    
    f = fopen([fname_ser,'coop',num2str(i),'.dat'],'r');
    A = mean(fscanf(f,'%f',[5,Inf]),2);
    fclose(f);
    timing_coop_ser.central = timing_coop_ser.central + A(1);
    timing_coop_ser.init = timing_coop_ser.init + ones(nmax,1)*[A(2), A(4)];
    timing_coop_ser.solve(i,:) = [A(3), A(5)];
    
    
    
    f = fopen([fname_par,'ncoop',num2str(i),'.dat'],'r');
    A = mean(fscanf(f,'%f',[5,Inf]),2);
    fclose(f);
    timing_ncoop_par.central = timing_ncoop_par.central + A(1);
    timing_ncoop_par.init = timing_ncoop_par.init + ones(nmax,1)*[A(2), A(4)];
    timing_ncoop_par.solve(i,:) = [A(3), A(5)];
    
    f = fopen([fname_ser,'ncoop',num2str(i),'.dat'],'r');
    A = mean(fscanf(f,'%f',[5,Inf]),2);
    fclose(f);
    timing_ncoop_ser.central = timing_ncoop_ser.central + A(1);
    timing_ncoop_ser.init = timing_ncoop_ser.init + ones(nmax,1)*[A(2), A(4)];
    timing_ncoop_ser.solve(i,:) = [A(3), A(5)];
    
end

timing_coop_par.central = timing_coop_par.central/nmax;
timing_coop_ser.central = timing_coop_ser.central/nmax;
timing_ncoop_par.central = timing_ncoop_par.central/nmax;
timing_ncoop_ser.central = timing_ncoop_ser.central/nmax;

timing_coop_par.init = timing_coop_par.init/nmax;
timing_coop_ser.init = timing_coop_ser.init/nmax;
timing_ncoop_par.init = timing_ncoop_par.init/nmax;
timing_ncoop_ser.init = timing_ncoop_ser.init/nmax;

timing_coop_par.solve = timing_coop_par.solve/nmax;
timing_coop_ser.solve = timing_coop_ser.solve/nmax;
timing_ncoop_par.solve = timing_ncoop_par.solve/nmax;
timing_ncoop_ser.solve = timing_ncoop_ser.solve/nmax;

timing_coop_par.total = timing_coop_par.central + ...
    mean(timing_coop_par.init,2) + mean(timing_coop_par.solve,2);

timing_coop_ser.total = timing_coop_ser.central + ...
    mean(timing_coop_ser.init,2) + mean(timing_coop_ser.solve,2);

timing_ncoop_par.total = timing_ncoop_par.central + ...
    mean(timing_ncoop_par.init,2) + mean(timing_ncoop_par.solve,2);

timing_ncoop_ser.total = timing_ncoop_ser.central + ...
    mean(timing_ncoop_ser.init,2) + mean(timing_ncoop_ser.solve,2);
