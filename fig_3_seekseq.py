import gzip
import concurrent.futures
import argparse
from datetime import datetime
import sys
import time
import Levenshtein as lv

class Logger:
    def __init__(self, log_file_name, stream=sys.stdout):
        self.log_file_name = log_file_name
        self.stream = stream

    def write(self, message):
        with open(self.log_file_name, 'a') as log_file:
            log_file.write(message)
        if self.stream:
            self.stream.write(message)

    def flush(self):
        if self.stream:
            self.stream.flush()

def merge_files(file_list, final_output):
    with open(final_output, 'w') as f_out:
        for fname in file_list:
            with open(fname) as f_in:
                for line in f_in:
                    f_out.write(line)
            os.remove(fname)  # Optionally remove the part file after merging

def process_block(r1_lines, r2_lines, const, scaffold, max_dist, index):
    len_const = len(const)
    len_sc = len(scaffold)
    matched_pairs = 0
    unique_events = set()

    fq1_output = open(f"part_fq1_{index}.txt", 'w')
    fq2_output = open(f"part_fq2_{index}.txt", 'w')

    while True:

        # Assume each block is already split into manageable parts
        r1_identifier = r1_lines[0]
        r1_sequence = r1_lines[2]
        r1_plus_line = r1_lines[3]
        r1_quality = r1_lines[4]

        print(r1_identifier)
        r2_identifier = r2_lines[0]
        r2_sequence = r2_lines[1]
        r2_plus_line = r2_lines[2]
        r2_quality = r2_lines[3]

        r1_match = False
        r2_match = False

        # Check R1
        for i in range(len(r1_sequence) - len_const + 1):
            if lv.distance(const, r1_sequence[i:i+len_const]) <= max_dist:
                r1_start_index = max(0, i-25)
                r1_end_index = min(len(r1_sequence), i + len_const + 10)
                event_sequence = r1_sequence[i+len_const:r1_end_index]
                r1_extracted_seq = r1_sequence[r1_start_index:i] + r1_sequence[i+len_const:r1_end_index]
                r1_extracted_quality = r1_quality[r1_start_index:i] + r1_quality[i+len_const:r1_end_index]
                if len(r1_extracted_seq) == 35:
                    r1_match = True
                break

        # Check R2
        for i in range(len(r2_sequence) - len_sc + 1, 0, -1): # Search from end would be faster
            if lv.distance(scaffold, r2_sequence[i:i+len_sc]) <= 1 and i-21 >= 0:
                r2_extracted_seq = r2_sequence[i-21:i]
                r2_extracted_quality = r2_quality[i-21:i]
                r2_match = True
                break
        
        if r1_match and r2_match:
            # check length
            # if len(r1_extracted_seq) == 35 and len(r2_extracted_seq) == 21: # redundant
                is_new_event = True
                for existing_event in unique_events:
                    if lv.distance(event_sequence, existing_event) <= 1:
                        is_new_event = False
                        break
                if is_new_event:
                    unique_events.add(event_sequence)
                    current_count += 1

                matched_pairs += 1
                unique_events.add(r1_identifier)  # Example of recording an event
                
                fq1_output.write(f"{r1_identifier}\n{r1_extracted_seq}\n{r1_plus_line}\n{r1_extracted_quality}\n")
                fq2_output.write(f"{r2_identifier}\n{r2_extracted_seq}\n{r2_plus_line}\n{r2_extracted_quality}\n")
        
        try:
            r1_lines = r1_lines[4:]
            r2_lines = r2_lines[4:]
        except IndexError:
            break
    return matched_pairs, len(unique_events)

def extract_paired_sequences(r1_input_file, fq1_output_file, r2_input_file, fq2_output_file, const, scaffold, max_dist=1, num_threads=4):
    
    log_file_name = f"seek_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
    logger = Logger(log_file_name)

    # Open files
    with gzip.open(r1_input_file, 'rt') as r1_file, open(fq1_output_file, 'wt') as fq1_output, \
         gzip.open(r2_input_file, 'rt') as r2_file, open(fq2_output_file, 'wt') as fq2_output:

        # Read all lines (for demonstration; in practice, read in chunks)
        r1_lines = [line.strip() for line in r1_file.readlines()]
        r2_lines = [line.strip() for line in r2_file.readlines()]

        # Divide work into chunks for each thread
        chunk_size = len(r1_lines) // 4 // num_threads
        futures = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:
            for i in range(num_threads):
                start_index = i * 4 * chunk_size
                end_index = start_index + 4 * chunk_size if i != num_threads - 1 else len(r1_lines)
                futures.append(executor.submit(process_block, r1_lines[start_index:end_index], r2_lines[start_index:end_index], const, scaffold, max_dist, i))

        # Collect results
        total_matched = 0
        total_unique = 0
        for future in concurrent.futures.as_completed(futures):
            matched, unique = future.result()
            total_matched += matched
            total_unique += unique

        # Merge part files
        for i in range(2):  # Two output files per thread
            part_files = [fname for j, fname in enumerate(output_files) if j % 2 == i]
            #merge_files(part_files, f"final_output_{i + 1}.txt")

    #logger.write(f"Finished processing. Total matched pairs: {total_matched}\n")
    #logger.write(f"Total unique events: {total_unique}\n")

# Command-line interface
if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="SeekSeq")
    parser.add_argument("-i", "--input", type=str, required=True, help="Input fastq File: <in1>,<in2>")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output fastq File: <out1>,<out2>")
    parser.add_argument("-c", "--constant", type=str, help="Constant Sequence in read1", default="TTGTCTTCCTAAGAC")
    parser.add_argument("-s", "--scaffold", type=str, help="Scaffold Sequence in read2", default="GTTTTAGA")
    parser.add_argument("-e", "--maxerror", type=int, help="Max wrong bases to match constant and scaffold", default=1)
    parser.add_argument("-@", "--threads", type=int, help="Max number of threads used", default=4)
    args = parser.parse_args()

    input_1, input_2 = args.input.split(',')
    output_1, output_2 = args.output.split(',')

    start_time = time.time()
    extract_paired_sequences(input_1, output_1, input_2, output_2, args.constant, args.scaffold, args.maxerror, args.threads)
    end_time = time.time()
    print("Total execution time: {:.2f} seconds".format(end_time - start_time))
