function S = calc_shannon(relative_abundance_vals)

    S = -1*sum(relative_abundance_vals.*log(relative_abundance_vals));

end