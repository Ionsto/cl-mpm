mcst = [9.92969028487503 0.304957091385821 0.1 0.1 2.92907060239431]
strain_in = [-15.5913110798307 9.76195134112692 2.59893240801835 19.0395407623297 4.09318197180819 -15.8508955283125]
eps = [strain_in(1),strain_in(2),strain_in(3),strain_in(6),strain_in(4),strain_in(5)]
strain_out = DPconst(eps,mcst);
swiz = [1 2 3 5 6 4]
strain_out(swiz)