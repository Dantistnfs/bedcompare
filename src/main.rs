use std::fs::File;
use std::io::{self, BufRead};
use clap::{Arg, Command};
use std::collections::{HashMap, HashSet};
use serde_json;

static EMPTY_BED_RECORDS: Vec<BedRecord> = vec![];

#[derive(Debug)]
struct BedRecord {
    start: u64,
    end: u64,
    gene: String,
}

fn read_bed_file(filename: &str) -> Result<HashMap<String,Vec<BedRecord>>, io::Error> {
    let file = File::open(filename).expect(&format!("Unable to open file: {}\n", filename));
    let reader = io::BufReader::new(file);

    let mut records = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split_whitespace().collect();

        if fields.len() >= 4 {
            let chromosome = fields[0].to_string();
            let record = BedRecord {
                start: fields[1].parse().expect(&format!("Malformed bed record in {}, line: {}", filename, line)),
                end: fields[2].parse().expect(&format!("Malformed bed record in {}, line: {}", filename, line)),
                gene: fields[3].to_string(),
            };

            let r = records.entry(chromosome).or_insert(Vec::new());
            r.push(record);
        }
    }

    Ok(records)
}

fn compare_genomic_files(file1: &str, file2: &str, json_output: bool) -> Result<(), io::Error> {
    let bed_records1 = read_bed_file(file1)?;
    let bed_records2 = read_bed_file(file2)?;

    let mut genome_length = 0;
    let mut total_bed_positions : u64 = 0;
    let mut same_gene_bp = 0;
    let mut same_regions = 0;
    let mut missing_regions = 0;
    let mut missing_chromosomes = 0;
    let mut number_of_regions = 0;
    let mut missing_chromosomes_file_1 = String::new();
    let mut missing_chromosomes_file_2 = String::new();


    let mut chromosomes : Vec<_> = bed_records1.keys().into_iter().collect();
    chromosomes.append(&mut bed_records2.keys().into_iter().collect());

    let chromosomes = chromosomes.into_iter().collect::<HashSet<_>>();

    for chromosome in chromosomes.iter() {
        let records1 = bed_records1.get(*chromosome).unwrap_or(&EMPTY_BED_RECORDS);
        let records2 = bed_records2.get(*chromosome).unwrap_or(&EMPTY_BED_RECORDS);

        let mut chromosome_length = 0;
        if let Some(rec) = records1.last() { chromosome_length = chromosome_length.max(rec.end); } 
        if let Some(rec) = records2.last() { chromosome_length = chromosome_length.max(rec.end); }
        genome_length += chromosome_length;

        if records1.len() == 0 {
            missing_chromosomes += 1;
            missing_chromosomes_file_1.push_str(chromosome);
            missing_chromosomes_file_1.push(',');
            let nregions = records2.iter().map(|rec| &rec.gene).collect::<HashSet<_>>().len();
            missing_regions += nregions;
            number_of_regions += nregions;
            total_bed_positions += records2.iter().map(|rec| rec.end.abs_diff(rec.start)).sum::<u64>();
            continue;
        }

        if records2.len() == 0 {
            missing_chromosomes += 1;
            missing_chromosomes_file_2.push_str(chromosome);
            missing_chromosomes_file_2.push(',');
            let nregions = records1.iter().map(|rec| &rec.gene).collect::<HashSet<_>>().len();
            missing_regions += nregions;
            number_of_regions += nregions;
            total_bed_positions += records1.iter().map(|rec| rec.end.abs_diff(rec.start)).sum::<u64>();
            continue;
        }

        let mut current_position = 0;
        let mut pointer1 = 0;
        let mut pointer2 = 0;
        // use pointer1 and pointer2 to move through records1 and records2, and calculate same and
        // different regions
        loop {
            let record1 = &records1[pointer1];
            let record2 = &records2[pointer2];
            //println!("lol kek");
            //dbg!(current_position);
            //dbg!(record1);
            //dbg!(record2);

            // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            //   record1.end   record2.end    current_position
            if record1.end <= current_position {
                if pointer1 < records1.len() - 1{
                    pointer1 += 1;
                    continue;
                }
            }
            // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            //   record2.end   record1.end    current_position
            if record2.end <= current_position {
                if pointer2 < records2.len() - 1{
                    pointer2 += 1;
                    continue;
                }
            }
            // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            //   record2.end   record1.end
            //                 current_position ------^
            //                    OR
            // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            //   record1.end   record2.end
            //                 current_position ------^
            if (record1.end <= current_position && record2.end == current_position) ||
               (record2.end <= current_position && record1.end == current_position) {
                // Try to exit loop updating records
                if (pointer1 == records1.len() - 1) && (pointer2 == records2.len() - 1) {
                    // Exit looping
                    break;
                }
                if pointer1 < records1.len() - 1{
                    pointer1 += 1;
                }
                if pointer2 < records2.len() - 1{
                    pointer2 += 1;
                }

                continue;
            }
            // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            //   current_position   record1.start   record2.start
            if current_position < record1.start.min(record2.start) {
                current_position = record1.start.min(record2.start);
            }
            // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            //  record1.start              record2.start
            //  current_position ----------^
            else if (current_position >= record1.start && current_position < record1.end) && current_position < record2.start {
                total_bed_positions += record2.start - current_position;
                current_position = record2.start;
            }
            // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            //  record2.start              record1.start
            //  current_position ----------^
            else if (current_position >= record2.start && current_position < record2.end) && current_position < record1.start {
                total_bed_positions += record1.start - current_position;
                current_position = record1.start;
            }
            // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            //  record1.start record2.start   record1.end record2.end
            //                current_position^
            else if (current_position >= record1.start && current_position >= record2.start) &&
                    (current_position < record1.end && current_position < record2.end) {
                    let smaller_end = record1.end.min(record2.end);
                    total_bed_positions += smaller_end - current_position;
                    number_of_regions += 1;

                    if record1.gene == record2.gene {
                        same_regions += 1;
                        same_gene_bp += smaller_end - current_position;
                    }
                    current_position = smaller_end;
            }
            // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            //  record1.end              record2.end
            //  current_position ----^
            else if current_position >= record1.end && current_position < record2.end && current_position >= record1.start {
                total_bed_positions += record2.end - current_position;
                current_position = record2.end;
            }
            // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            //  record2.end              record1.end
            //  current_position ----^
            else if current_position >= record2.end && current_position < record1.end && current_position >= record1.start {
                total_bed_positions += record1.end - current_position;
                current_position = record1.end;
            }
            // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            //  record1.start record2.start record1.end record2.end
            //                                  current_position-----^
            else if (current_position >= record1.start && current_position >= record2.start) &&
                    (current_position > record1.end && current_position > record2.end) {
                current_position = record1.start.min(record2.start);
            }
            // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            //  record1.start record1.end    record2.start record2.end
            //       current_position-----^
            else if (current_position >= record1.start && current_position >= record1.end) &&
                    (current_position <= record2.start && current_position <= record2.end){
                current_position = record2.start;
            }
            // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            //  record2.start record2.end    record1.start record1.end
            //       current_position-----^
            else if (current_position >= record2.start && current_position >= record2.end) &&
                    (current_position <= record1.start && current_position <= record1.end){
                current_position = record1.start;
            }

            continue;
        }

    }

    let same_region_bp_percent = same_gene_bp as f64/total_bed_positions as f64;
    
    if json_output {
        let mut output = HashMap::new();
        output.insert("genome_length", format!("{}", genome_length as f64));
        output.insert("bed_regions_length", format!("{}", total_bed_positions as f64));
        output.insert("bed_number_of_regions", format!("{}", number_of_regions as f64));
        output.insert("same_regions", format!("{}", same_regions as f64));
        output.insert("same_regions_length", format!("{}", same_gene_bp as f64));
        output.insert("same_regions_length_percent", format!("{}", same_region_bp_percent as f64));
        output.insert("missing_chromosomes", format!("{}", missing_chromosomes as f64));
        output.insert("missing_regions", format!("{}", missing_regions as f64));
        output.insert("missing_chromosomes_file_1", missing_chromosomes_file_1);
        output.insert("missing_chromosomes_file_2", missing_chromosomes_file_2);
        println!("{}", serde_json::to_string(&output).unwrap());
    }
    else {
    println!("Assumed genome length: {}", genome_length);
    println!("Length of geneome in bed files regions: {}", total_bed_positions);
    println!("Total number of regions: {}", number_of_regions);

    println!("Number of same regions: {}", same_regions);
    println!("Length of same regions: {}", same_gene_bp);
    println!("Percent of same bp regions: {}%", same_region_bp_percent * 100.0);

    println!("Chromosomes missing from one of the files: {}", missing_chromosomes);
    println!("Number of regions missing from files: {}", missing_regions);
    }
                                        
    Ok(())
}


fn main() {
    let matches = Command::new("bedcompare")
        .version("0.1")
        .author("Stanislav Zubenko")
        .about("Compares two genomic .bed files on base pair level")
        .arg_required_else_help(true)
        .arg(
            Arg::new("file1")
                .help("Path to the first .bed file")
                .short('1')
                .long("bed1")
                .required(true)
        )
        .arg(
            Arg::new("file2")
                .help("Path to the second .bed file")
                .short('2')
                .long("bed2")
                .required(true)
        )
        .arg(
            Arg::new("json")
                .num_args(0)
                .help("Output in json format")
                .short('j')
                .long("json")
        )
        .get_matches();

    let file1 = matches.get_one::<String>("file1").unwrap();
    let file2 = matches.get_one::<String>("file2").unwrap();
    let json_output = matches.get_one::<bool>("json").unwrap();
    if let Err(err) = compare_genomic_files(file1, file2, *json_output) {
        eprintln!("Error: {}", err);
    }
}
