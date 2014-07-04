ls Gen*/best_permutation_distribution_hyp.png | sort -k2 -tn -n > gif_loc_hyp.txt
ls Gen*/best_permutation_ratios_0.1Kdiv.png | sort -k2 -tn -n > gif_loc_rat.txt # add command line argument

convert @gif_loc_hyp.txt images_hyp.gif
convert @gif_loc_rat.txt images_ratios.gif

ls Gen*/best_permutation_distribution_hm.png | sort -k2 -tn -n > gif_loc_hm.txt
convert @gif_loc_hm.txt images_hm.gif