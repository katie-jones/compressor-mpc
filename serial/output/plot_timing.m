N = 9;
n_controllers = 2;
n_iterations = 250/0.05;

t_init = zeros(N,n_controllers);
t_solve = t_init;

for i=1:N
    f = fopen(['ncoop_timing',num2str(i),'.dat'],'r');
    A = fscanf(f,'%g');
    fclose(f);
    t_init(i,:) = A([1,3])/n_iterations/1.0e6;
    t_solve(i,:) = A([2,4])/n_iterations/1.0e6;
end