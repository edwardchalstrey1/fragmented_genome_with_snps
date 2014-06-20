ls Gen*/best_permutation_distribution_hyp.png | sort -k2 -tn -n > gif_loc_hyp.txt
ls Gen*/best_permutation_ratios_10Kdiv.png | sort -k2 -tn -n > gif_loc_rat.txt

convert @gif_loc_hyp.txt images_hyp.gif
convert @gif_loc_rat.txt images_ratios.gif