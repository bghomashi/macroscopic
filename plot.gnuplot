reset

set logscale y; 
set yrange [1e-15:1];
plot [0:80]  "he_data_lerp_1e0.txt" u ($1/0.057):($2**2 + $3**2) w l lw 2, \
             "he_data_lerp_1e4.txt" u ($1/0.057):($2**2 + $3**2) w l lw 2, \
             "he_data_lerp_1e5.txt" u ($1/0.057):($2**2 + $3**2) w l lw 2, \
             "he_data_lerp_1e6.txt" u ($1/0.057):($2**2 + $3**2) w l lw 2

             
#plot [0:23]  "joel_data_lerp_2_1e0.txt" u ($1/0.057):($2**2 + $3**2) w l lw 2, \
#             "joel_data_lerp_2_1e4.txt" u ($1/0.057):($2**2 + $3**2) w l lw 2