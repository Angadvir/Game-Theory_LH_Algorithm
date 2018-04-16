global en totC counter count totCo;
totC = 0; totCo =0;counter =0;k =[2 3 4 5 6 7 8 9 10]; 
clc;
%PLEASE ALLOW CODE TO RUN FOR APPROX. 145 SECONDS.
%IF OUTPUT SHOWS, 'DEGENERACY', PLEASE RUN ONCE AGAIN.
for en=2:10 
    totCo = 0;
    for i=1:206
     LH2_part1();
     len = length(LH2_part3__2());
     totCo= totCo+count;
     totC(en-1)= totCo/206;
    end  
   fprintf(1,'Thanks for waiting up until iteration: %d\n',en);   
end 
scatter(k, totC);
xlabel('Size of "n" matrix');
ylabel('Number of pivots to find a solution/1000 iterations');