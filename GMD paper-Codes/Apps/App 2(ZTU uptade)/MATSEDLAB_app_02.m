rt_time = clock;                                    %display starting time
fprintf('\nThe starting time is %.0f/%.0f/%.0f %.0f:%.0f:%.0f\n', rt_time(1),...
    rt_time(2), rt_time(3), rt_time(4), rt_time(5), rt_time(6));

MATSEDLAB;

global n_segment;
n_segment = 2;
t = 100;
for i1 = 1:100
    disp(t);
    MATSEDLAB_ox;
    MATSEDLAB_anox;
    t = t+1;
end

global SimValues para;
save('Result.mat','SimValues','para');        %save the result in Result.mat

rt_time = clock;                                    %display ending time
fprintf('\nThe ending time is %.0f/%.0f/%.0f %.0f:%.0f:%.0f\n', rt_time(1),...
    rt_time(2), rt_time(3), rt_time(4), rt_time(5), rt_time(6));