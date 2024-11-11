#collect QUAST results for all genomes in one file
#usage: bash quast_results.sh

head -1 quast_results/GCA_000007465_quast/transposed_report.tsv > quast_results.tsv

for report_file in quast_results/*/transposed_report.tsv
        do
        head -2 $report_file | tail -1 >> quast_results.tsv
	done
