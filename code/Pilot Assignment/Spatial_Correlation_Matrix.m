function  R_corr = Spatial_Correlation_Matrix(A, NUM_USERS,NUM_APs,AP_ANTENNAS)
R_corr          = zeros(AP_ANTENNAS,AP_ANTENNAS,NUM_APs,NUM_USERS) ;
for k = 1:NUM_USERS
    for l = find(A(:,k))'
        for i = 1:AP_ANTENNAS
            for j = 1:AP_ANTENNAS
                
                R_corr(i,j,l,k) = 1;
                
            end
        end
    end
end

end