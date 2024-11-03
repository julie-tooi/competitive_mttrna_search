# the purpose of this script: add information about codon in annotation
# search for the anticodon loop (by structure)
# add the info about codon (can be NA)
# add the information about one stranded loop element for next group step

import re
from collections import defaultdict
import argparse
from operator import itemgetter


def translate(dna_seq):
    translate_table = {'A': 'U', 'T': 'A', 'C': 'G', 'G': 'C'}

    rna_seq = []

    for base in dna_seq:
        rna_seq.append(translate_table[base])

    return ''.join(rna_seq)


def collect_trna_anticodons(anticodon_collection):
    # work with the NCBI codes
    trna_anticodons = defaultdict(set)

    # collect codons in trna types
    with open(anticodon_collection, 'r') as file_handler:
        for line in file_handler:
            amino_acids = line.strip().split(" = ")[-1]
            file_handler.readline()  # starts
            base_1 = file_handler.readline().strip().split(" = ")[-1]
            base_2 = file_handler.readline().strip().split(" = ")[-1]
            base_3 = file_handler.readline().strip().split(" = ")[-1]

        for amino_acid_pos in range(len(amino_acids)):
            trna_gene = f"trn{amino_acids[amino_acid_pos]}"
            codon = f"{base_1[amino_acid_pos]}" \
                    f"{base_2[amino_acid_pos]}" \
                    f"{base_3[amino_acid_pos]}"

            anticodon = translate(codon)[::-1]

            trna_anticodons[trna_gene].add(anticodon)

    return trna_anticodons


def update_structure_and_seq(sequence, structure):
    allowed_bases = ["A", "U", "G", "C"]
    new_structure = []
    new_seq = []
    for base, base_str in zip(sequence, structure):
        if base.upper() in allowed_bases:
            new_structure.append(base_str)
            new_seq.append(base)

    new_structure = "".join(new_structure)
    new_seq = "".join(new_seq).upper()

    return new_structure, new_seq


def find_loops(structure):
    inner_loops_search = re.compile(r"<([.,~_-]*)>")

    search_result_loops = inner_loops_search.finditer(structure)

    return list(search_result_loops)


def define_ac_loop(loops, shape):
    new_loops = []
    for loop in loops:
        new_loops.append(loop.group(1))
    number_of_loops = len(new_loops)

    anticodon_loop = None
    defined_shape = "unclear"

    if (shape == "cloverleaf") and (number_of_loops == 3):
        if len(new_loops[1]) > 6:
            anticodon_loop = loops[1]
            defined_shape = "D_loop_anticodon_T_loop"

    elif (shape == "one-armless") and (number_of_loops == 2):
        if (len(new_loops[0]) > 6) and (len(new_loops[1]) > 6):
            anticodon_loop = [loops[0], loops[1]]
            defined_shape = ["anticodon_T_loop", "D_loop_anticodon"]
        elif len(new_loops[0]) > 6:
            anticodon_loop = loops[0]
            defined_shape = "anticodon_T_loop"
        elif len(new_loops[1]) > 6:
            anticodon_loop = loops[1]
            defined_shape = "D_loop_anticodon"

    elif (shape == "two-armless") and (number_of_loops == 1):
        if len(new_loops[0]) > 6:
            anticodon_loop = loops[0]
            defined_shape = "anticodon"

    return anticodon_loop, defined_shape


def search_codon(codon_loop, sequence, trn_type, anticodon_collection):
    codon_loop_start, codon_loop_end = codon_loop.start(1), \
                                       codon_loop.end(1)

    loop_sequence = sequence[codon_loop_start:codon_loop_end + 1]

    anticodon_candidates = []

    if len(trn_type) == 5:
        trn_type = trn_type[:-1]

    for anticodon in anticodon_collection[trn_type]:
        anticodon_search_results = list(re.finditer(anticodon, loop_sequence))
        if anticodon_search_results:
            for anticodon in anticodon_search_results:
                if anticodon.start(0) >= 1 and \
                        anticodon.end(0) <= (len(loop_sequence) - 1):
                    anticodon_candidates.append(anticodon)

    return anticodon_candidates


def define_groups(strand_group):
    sorted_strand_groups = sorted(strand_group, key=itemgetter(0))

    one_group = []
    groups_collection = []

    for group in sorted_strand_groups:
        if not one_group:
            one_group.append(group)
        else:
            # check coordinate intersections:
            previous_group_start = int(one_group[-1][0])
            previous_group_end = int(one_group[-1][1])

            group_start = group[0]
            group_end = group[1]

            if (previous_group_start <= group_start <= previous_group_end) or \
                    (previous_group_start <= group_end <= previous_group_end):
                one_group.append(group)
            else:
                groups_collection.append(one_group)
                one_group = [group]

    if one_group:
        groups_collection.append(one_group)

    return groups_collection


def find_top_hit(one_group):
    sorted_by_bit_score = sorted(one_group, key=lambda start: float(start[2][2]))

    top_hit = sorted_by_bit_score[-1]

    return top_hit


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, help='Path to collected annotation', required=True)
    parser.add_argument('-o', '--output', type=str, help='Path to updated and filtered annotation', required=True)
    parser.add_argument('-c', '--codon', action='store_true', help='Require codon table: strict codon filter')
    parser.add_argument('-t', '--trn', type=str, help='Path to trna-anticodon table (NCBI format)', required=True)
    parser.add_argument('-f', '--filter_threshold', type=int, help='Bit-score filter: set the filter threshold',
                        required=True)

    args = parser.parse_args()

    path_to_annnotation = args.input
    path_to_filtered_output = args.output
    path_to_anticodon_collection = args.trn
    strict_codon_filter = args.codon
    filter_threshold = args.filter_threshold

    trna_anticodon_collection = collect_trna_anticodons(path_to_anticodon_collection)

    with open(path_to_filtered_output, 'w') as output:
        with open(path_to_annnotation, 'r') as collection_handler:

            header = collection_handler.readline().strip().split("\t")
            header.extend([
                "anticodon", "ac_loop_start",
                "ac_loop_end", "shape"])
            header = "\t".join(header)

            output.write(f"{header}\n")

            # collect annotation info
            for line in collection_handler:
                line = line.strip().split("\t")

                seq_start, seq_end = int(line[3]), int(line[4])
                initial_structure = line[6]
                initial_shape = line[8]
                trna_gene = line[9]
                target_sequence = line[10]

                updated_structure, complete_seq = update_structure_and_seq(target_sequence,
                                                                           initial_structure)

                inner_loops = find_loops(updated_structure)

                # ac loop before checks:search by structure and length
                anticodon_loop, trna_shape = define_ac_loop(inner_loops,
                                                            initial_shape)

                if anticodon_loop:
                    if type(anticodon_loop) == list:
                        anticodon_coordinates_loop_1 = search_codon(anticodon_loop[0],
                                                                    complete_seq,
                                                                    trna_gene,
                                                                    trna_anticodon_collection)
                        anticodon_coordinates_loop_2 = search_codon(anticodon_loop[1],
                                                                    complete_seq,
                                                                    trna_gene,
                                                                    trna_anticodon_collection)

                        # if we see anticodon in both loops:
                        # we cannot say anything about shape but save both entry
                        if anticodon_coordinates_loop_1 and anticodon_coordinates_loop_2:
                            trna_shape = "unclear"
                            codon_search_result = f"{anticodon_coordinates_loop_1[0].group(0)}"
                            loop_start = anticodon_loop[0].start(1) + seq_start
                            loop_end = anticodon_loop[0].end(1) + seq_start
                            # first entry
                            line_first_entry = line
                            line_first_entry.extend([codon_search_result,
                                                     str(loop_start), str(loop_end),
                                                     trna_shape])
                            line_first_entry = "\t".join(line_first_entry)
                            output.write(f"{line_first_entry}\n")

                            # prepare data for second entry
                            codon_search_result = f"{anticodon_coordinates_loop_2[0].group(0)}"
                            loop_start = anticodon_loop[1].start(1) + seq_start
                            loop_end = anticodon_loop[1].end(1) + seq_start

                        # one loop anticodon case
                        elif anticodon_coordinates_loop_1:
                            trna_shape = "anticodon_T_loop"
                            codon_search_result = anticodon_coordinates_loop_1[0].group(0)
                            loop_start = anticodon_loop[0].start(1) + seq_start
                            loop_end = anticodon_loop[0].end(1) + seq_start

                        elif anticodon_coordinates_loop_2:
                            trna_shape = "D_loop_anticodon"
                            codon_search_result = anticodon_coordinates_loop_2[0].group(0)
                            loop_start = anticodon_loop[1].start(1) + seq_start
                            loop_end = anticodon_loop[1].end(1) + seq_start

                        else:
                            trna_shape = "unclear"
                            codon_search_result = "N/A"
                            loop_start = anticodon_loop[0].start(1) + seq_start
                            loop_end = anticodon_loop[1].end(1) + seq_start

                    else:
                        anticodon_coordinates = search_codon(anticodon_loop,
                                                             complete_seq,
                                                             trna_gene,
                                                             trna_anticodon_collection)

                        if anticodon_coordinates:
                            codon_search_result = anticodon_coordinates[0].group(0)
                            loop_start = anticodon_loop.start(1) + seq_start
                            loop_end = anticodon_loop.end(1) + seq_start
                        else:
                            codon_search_result = "N/A"
                            loop_start = anticodon_loop.start(1) + seq_start
                            loop_end = anticodon_loop.end(1) + seq_start

                else:
                    codon_search_result = "N/A"
                    loop_start, loop_end = "N/A", "N/A"

                line.extend([codon_search_result,
                             str(loop_start), str(loop_end),
                             trna_shape])
                line = "\t".join(line)
                output.write(f"{line}\n")

    # define locus group
    pos_strand_group = []
    neg_strand_group = []

    with open(path_to_filtered_output, 'r') as prefiltered_handler:
        header = prefiltered_handler.readline().strip().split("\t")  # header
        header.extend((
            "group", "top_hit"
        ))

        header = "\t".join(header)

        for line in prefiltered_handler:
            line = line.strip().split("\t")
            anticodon, ac_loop_start, ac_loop_end = line[11], line[12], line[13]
            score = float(line[2])

            strand = line[5]

            if score > filter_threshold and ac_loop_start != "N/A":
                if strict_codon_filter:
                    if anticodon != "N/A":
                        ac_loop_start, ac_loop_end = int(ac_loop_start), \
                                                     int(ac_loop_end)
                        if strand == "+":
                            pos_strand_group.append((
                                ac_loop_start, ac_loop_end, line
                            ))
                        else:
                            neg_strand_group.append((
                                ac_loop_start, ac_loop_end, line
                            ))
                else:
                    ac_loop_start, ac_loop_end = int(ac_loop_start), \
                                                 int(ac_loop_end)
                    if strand == "+":
                        pos_strand_group.append((
                            ac_loop_start, ac_loop_end, line
                        ))
                    else:
                        neg_strand_group.append((
                            ac_loop_start, ac_loop_end, line
                        ))

        pos_groups = define_groups(pos_strand_group)
        neg_groups = define_groups(neg_strand_group)

        with open(f"{path_to_filtered_output}_locus", 'w') as collection_handler:
            collection_handler.write(f"{header}\n")

            for index_group, group in enumerate(pos_groups):

                group_name = f"{index_group}+"
                top_hit = find_top_hit(group)

                for candidate in group:
                    if candidate == top_hit:
                        new_line = candidate[2]
                        new_line.extend((
                            group_name, "True"
                        ))
                    else:
                        new_line = candidate[2]
                        new_line.extend((
                            group_name, "False"
                        ))

                    new_line = "\t".join(new_line)
                    collection_handler.write(f"{new_line}\n")

            for index_group, group in enumerate(neg_groups):

                group_name = f"{index_group}-"
                top_hit = find_top_hit(group)

                for candidate in group:
                    if candidate == top_hit:
                        new_line = candidate[2]
                        new_line.extend((
                            group_name, "True"
                        ))
                    else:
                        new_line = candidate[2]
                        new_line.extend((
                            group_name, "False"
                        ))

                    new_line = "\t".join(new_line)
                    collection_handler.write(f"{new_line}\n")

    print("Done")
