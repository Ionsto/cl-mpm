size = 1000
all_data = zeros(size,6+6+5);
for i = 1 : size
    %E,nu,phi,psi,c
    phi = rand
    if phi == 0
        phi = 0.1
    end
    mcst = [rand * 10,rand*0.4999,.1,0.1,rand*10];
    strain_in = (rand(6,1) * 50) - 25;
    strain_out = DPconst(strain_in,mcst);
    all_data(i,:) = [strain_in',mcst,strain_out']';
end
MTR_titles = {'in_e_xx' 'in_e_yy' 'in_e_zz' 'in_e_xy' 'in_e_yz' 'in_e_zx' 'E' 'nu' 'phi' 'psi' 'c' 'out_e_xx' 'out_e_yy' 'out_e_zz' 'out_e_xy' 'out_e_yz' 'out_e_zx'};
all_data = array2table(all_data,'VariableNames',MTR_titles);
writetable(all_data,"dp_test_data.csv");
