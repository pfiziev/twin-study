rm CGmap/log
python annotate_regions.py
python annotate_sites.py
python calc_fragment_methylation_levels.py
python compute_average_methylation.py

# After normalization
python calc_average_methylation_per_region_on_fragment_basis.py
python calc_correlation.py
python find_fragments_for_age_prediction.py
python regression.py

