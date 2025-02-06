
function printBeautyMatrix(matr1)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[r_count, c_count]=size(matr1);
matr_name=inputname(1);
fprintf('%s = left (matrix {',matr_name);
for r=1:r_count
    for c=1:c_count
        if (c~=c_count)
            fprintf("%.6f #", matr1(r,c));
        else
            fprintf("%.6f", matr1(r,c)); 
        end
    end
    if (r~=r_count)
        fprintf("## ");
    end
end
fprintf('} right )\n');
end
