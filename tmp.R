study = trscr_lusc
targeted_genes = penda_superdown_deregulated
index_pathology = 2
up_str = 2500
dwn_str = 2500
nb_probe_min = 1


plot_selected_gene("GRK5",expr_data = trscr_lusc$data, meth_data = meth_lusc$data, probes_index = features$probes)
plot_selected_probes("GRK5",nbins = 6)