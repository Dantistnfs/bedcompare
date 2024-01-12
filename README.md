# Bedcompare

Tool for comparison of genomic bed files

It assumes two things about bed files:
- bed files are sorted
- regions in bed files are not intersecting

### Usage
```bash
> bedcompare -1 ./bed/variants_68_percent.bed -2 ./bed/variants_68_percent_1.bed --json
{"bed_regions_length":1899.0,"same_regions_length":1299.0,"genome_length":2000.0,"same_regions":3.0,"missing_chromosomes":0.0,"missing_regions":0.0,"same_regions_length_percent":0.684044233807267,"bed_number_of_regions":3.0}
```
