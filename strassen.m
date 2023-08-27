function [A,B,C] = strassen() 
% generates Strassen's factor matrices for 2x2x2

    A = [
      1  0  1  0  1 -1  0 
      0  1  0  0  0  1  0 
      0  0  0  0  1  0  1 
      1  1  0  1  0  0 -1
    ];
    
    B = [
      1  1  0 -1  0  1  0 
      0  0  0  1  0  0  1 
      0  0  1  0  0  1  0 
      1  0 -1  0  1  0  1     
    ];
    
    C = [
      1  0  0  1 -1  0  1 
      0  0  1  0  1  0  0 
      0  1  0  1  0  0  0 
      1 -1  1  0  0  1  0   
    ];

end