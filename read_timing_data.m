nmax = 9;
nruns = 5;

timing_coop_par = zeros(nmax,nruns);
timing_ncoop_par = timing_coop_par;
timing_centralized_par = timing_coop_par;
timing_coop_ser = timing_coop_par;
timing_ncoop_ser = timing_ncoop_par;
timing_centralized_ser = timing_coop_par;

fname_par = 'parallel/';
fname_ser = 'serial/';

%% centralized
for run=1:nruns
    folder_name = ['run',num2str(run),'/'];
f = fopen([fname_par, folder_name, 'centralized.dat'],'r');
A = mean(fscanf(f,'%f',[21,Inf]),2);
fclose(f);
timing_centralized_par(:,run) = A(end)*ones(nmax,1);

f = fopen([fname_ser, folder_name, 'centralized.dat'],'r');
A = mean(fscanf(f,'%f',[20,Inf]),2);
fclose(f);
timing_centralized_ser(:,run) = A(end)*ones(nmax,1);

%% cooperative/non-cooperative

for i=1:nmax
    f = fopen([fname_par, folder_name, 'coop',num2str(i),'.dat'],'r');
    A = mean(fscanf(f,'%f',[21,Inf]),2);
    fclose(f);
    timing_coop_par(i,run) = A(end)/2;
    
    f = fopen([fname_ser, folder_name, 'coop',num2str(i),'.dat'],'r');
    A = mean(fscanf(f,'%f',[20,Inf]),2);
    fclose(f);
    timing_coop_ser(i,run) = A(end)/2;
    
    
    f = fopen([fname_par, folder_name, 'ncoop',num2str(i),'.dat'],'r');
    A = mean(fscanf(f,'%f',[21,Inf]),2);
    fclose(f);
    timing_ncoop_par(i,run) = A(end)/2;
    
    f = fopen([fname_ser, folder_name, 'ncoop',num2str(i),'.dat'],'r');
    A = mean(fscanf(f,'%f',[20,Inf]),2);
    fclose(f);
    timing_ncoop_ser(i,run) = A(end)/2;
    
end

end

res.parallel.ncoop = mean(timing_ncoop_par,2);
res.parallel.coop = mean(timing_coop_par,2);
res.parallel.cent = mean(timing_centralized_par,2);

res.serial.ncoop = mean(timing_ncoop_ser,2);
res.serial.coop = mean(timing_coop_ser,2);
res.serial.cent = mean(timing_centralized_ser,2);
