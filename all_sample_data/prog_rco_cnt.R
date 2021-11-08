filename = "progeny_recom_count_data_for_distribution_comparison.csv"
df = read.csv( filename, header = TRUE )

print( head(df))

wilcox.test( x = df$SCx529L, y = df$SCxP6 )