%% Prime Vector
INITIAL_MATRIX = cell(6, 1);
for CUR = 4:6
    
    sz = length(Prime_Vector{CUR});
    
    RNG_Vector = zeros(sz, 1);
    ABS_Error_Vector = zeros(sz, 1);
    REL_Error_Vector = zeros(sz, 1);
    RND_Error_Vector = zeros(sz, 1);
    FCNVAL_Vector = zeros(sz, 1);
    NUMITR_Vector = zeros(sz, 1);
    LOC_Vector = zeros(sz, 1);
    
    for i = 1:sz
        LOC_Vector(i) = i;
        RNG_Vector(i) = Prime_Vector{CUR}{i}.num_val;
        ABS_Error_Vector(i) = Prime_Vector{CUR}{i}.errors.abs;
        REL_Error_Vector(i) = Prime_Vector{CUR}{i}.errors.rel;
        RND_Error_Vector(i) = Prime_Vector{CUR}{i}.errors.rnd;
        FCNVAL_Vector(i) = Prime_Vector{CUR}{i}.out_struct.FcnVal;
        NUMITR_Vector(i) = Prime_Vector{CUR}{i}.out_struct.NumIter;
    end

    INITIAL_MATRIX{CUR} = [LOC_Vector RNG_Vector NUMITR_Vector FCNVAL_Vector ABS_Error_Vector REL_Error_Vector RND_Error_Vector];       
end
clear CUR i LOC_Vector RNG_Vector NUMITR_Vector FCNVAL_Vector ABS_Error_Vector REL_Error_Vector RND_Error_Vector;
%% SP_Data
RSP_MATRIX = cell(6, 1);

for CUR = 4:6
    sz = length(Prime_Vector{CUR});
    ABS_Error_Matrix = zeros(sz, 5, 7);
    REL_Error_Matrix = zeros(sz, 5, 7);
    RND_Error_Matrix = zeros(sz, 5, 7);
        
    for i = 1:7
        for j = 1:sz
            for k = 1:5
                ABS_Error_Matrix(j, k, i) = SP_Data{CUR}{j, k, i}.rsp_errors.abs;
                REL_Error_Matrix(j, k, i) = SP_Data{CUR}{j, k, i}.rsp_errors.rel;
                RND_Error_Matrix(j, k, i) = SP_Data{CUR}{j, k, i}.rsp_errors.rnd;
            end
        end

    end
    RSP_MATRIX{CUR} = [ABS_Error_Matrix REL_Error_Matrix RND_Error_Matrix];
end

clear i j k CUR sz ABS_Error_Matrix REL_Error_Matrix RND_Error_Matrix;
%% CP_Data

CP_MATRIX = cell(6, 1);

for CUR = 4:6
    sz = length(Prime_Vector{CUR});

    NUMITR_Matrix = zeros(sz, 1, 7);
    FCNVAL_Matrix = zeros(sz, 1, 7);
    ABS_Error_Matrix = zeros(sz, 5, 7);
    REL_Error_Matrix = zeros(sz, 5, 7);
    RND_Error_Matrix = zeros(sz, 5, 7);
    for i = 1:7
        for j = 1:sz
            for k = 1:5
                NUMITR_Matrix(j, k, i) = CP_Data{CUR}{j, k, i}.cp_output.NumIter;
                FCNVAL_Matrix(j, k, i) = CP_Data{CUR}{j, k, i}.cp_output.FcnVal;
                ABS_Error_Matrix(j, k, i) = CP_Data{CUR}{j, k, i}.cp_errors.abs;
                REL_Error_Matrix(j, k, i) = CP_Data{CUR}{j, k, i}.cp_errors.rel;
                RND_Error_Matrix(j, k, i) = CP_Data{CUR}{j, k, i}.cp_errors.rnd;
            end
        end
    end
    CP_MATRIX{CUR} = [NUMITR_Matrix FCNVAL_Matrix ABS_Error_Matrix REL_Error_Matrix RND_Error_Matrix];
end

clear i j k NUMITR_Matrix FCNVAL_Matrix CUR sz ABS_Error_Matrix REL_Error_Matrix RND_Error_Matrix;