function output = singularity_check(matrix)
[x,y] = size(matrix);
for i = 1:x
    for j = 1:y
        if matrix(i,j) == 0
            matrix(i,j) = 0.01;
        end
    end
end
output = matrix;