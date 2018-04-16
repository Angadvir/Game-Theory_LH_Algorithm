global en totC counter ;
totC = 0; counter =0; anss = 0;k =[2 3 4 5 6 7 8 9 10];
%PLEASE ALLOW CODE TO RUN FOR APPROX. 95 SECONDS.
for en=2:10 
    counter = 0
    for i=1:1000
     LH2_part1();
     len = length(LH2_part1());
     anss = LH2_part1();
     if anss == 1
      counter = counter + 1;
     end
     totC(en-1)= counter;
    end  
    
end 
scatter(k, totC);
xlabel('Size of "n" matrix');
ylabel('Number of Secondary Rays encountered/1000 runs');