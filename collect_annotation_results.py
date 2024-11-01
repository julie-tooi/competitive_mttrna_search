# collect annotation results from CM search

import argparse
import re
from pathlib import Path


def find_loops(structure):
    acceptor_stem_search = re.compile(r"(\(+)[<>.,~_-]*(\)+)")
    inner_loops_search = re.compile(r"(<+)[.,~_-]*(>+)")

    search_result_loops = inner_loops_search.findall(structure)
    search_result_stem = acceptor_stem_search.findall(structure)

    # check if acceptor stem annotated the same as others (<>)
    # sometimes it can happen with two-armless CM annotation
    if (not re.search("\(", structure)) and (not re.search("\)", structure)):
        if re.search("<", structure) and re.search(">", structure):
            search_result_loops = inner_loops_search.findall(structure)
            search_result_stem = None

    return search_result_loops, search_result_stem


def structure_description(loops, stem):
    if stem and len(loops) >= 3:
        description = "cloverleaf"
    elif stem and len(loops) == 2:
        description = "one-armless"
    elif stem and len(loops) == 1:
        description = "two-armless"
    else:
        description = "unclear"

    return description


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, help='Path to folder with the annotation', required=True)
    parser.add_argument('-o', '--output', type=str, help='Path to table with results', required=True)

    args = parser.parse_args()

    path_to_annotation_batch = Path(args.input)
    path_to_output = Path(args.output)

    files_batch = path_to_annotation_batch.glob("*_annotation")

    with open(path_to_output, 'w') as output:
        header = "\t".join(["mitogenome_id",
                            "eval", "score",
                            "start", "end", "strand",
                            "predicted_structure", "simple_structure",
                            "shape_description",
                            "cm_name", "trna_type",
                            "sequence_w_gaps"])
        output.write(f"{header}\n")

        for file in files_batch:

            with open(file, 'r') as file_handler:
                for line in file_handler:

                    # main parsing block
                    if line.startswith(">> "):
                        line_construct = []

                        parsed_line = line.strip().split()
                        candidate_id = parsed_line[1]

                        file_handler.readline() # header
                        file_handler.readline() # spacer

                        stats = file_handler.readline().strip().split()

                        evalue, score, start, end, strand = stats[2], stats[3], stats[9], stats[10], stats[11]

                        file_handler.readline() # empty line
                        file_handler.readline() # NC negative scoring non-canonical basepairs
                        structure = file_handler.readline().strip().split()[0] # CS predicted secondary structure of the target sequence

                        loops, stem = find_loops(structure)
                        structure_shape = structure_description(loops, stem)

                        updated_structure = structure.replace('<', '(').\
                            replace('>', ')').\
                            replace(',', '_').\
                            replace('.', '_').\
                            replace('-', '_')

                        model_info = file_handler.readline().strip().split()[0] # consensus of the query model, the highest scoring residue sequence is shown

                        file_handler.readline() # scoring data

                        sequence = file_handler.readline().strip()
                        target_sequence = sequence.split()[2]  # target sequence
                        target_sequence_start, target_sequence_end = sequence.split()[1], sequence.split()[-1]

                        file_handler.readline()  # empty line

                        # check if several lines in alignment
                        if (target_sequence_start != start) or (target_sequence_end != end):
                            file_handler.readline() # NC negative scoring non-canonical basepairs

                            additional_structure = file_handler.readline().strip().split()[0]
                            # CS predicted secondary structure of the target sequence
                            structure = structure + additional_structure

                            add_loops, add_stem = find_loops(structure)
                            structure_shape = structure_description(add_loops, add_stem)

                            updated_structure = structure.replace('<', '('). \
                                replace('>', ')'). \
                                replace(',', '_'). \
                                replace('.', '_')

                            file_handler.readline() # consensus of the query model, the highest scoring residue sequence is shown
                            file_handler.readline()  # scoring data

                            additional_sequence = file_handler.readline().strip()
                            target_sequence = target_sequence + additional_sequence.split()[2]  # full target sequence
                            target_sequence_start, target_sequence_end = additional_sequence.split()[1], \
                                                                         additional_sequence.split()[3]

                            if target_sequence[-1] == '[':
                                target_sequence = target_sequence + " " + additional_sequence.strip().split()[3]
                            if target_sequence[-1] == '[':
                                target_sequence = target_sequence + " " + additional_sequence.strip().split()[4]
                            if target_sequence[-1] == '[':
                                target_sequence = target_sequence + " " + additional_sequence.strip().split()[5]
                            if target_sequence[-1] == '[':
                                target_sequence = target_sequence + " " + additional_sequence.strip().split()[6]

                            file_handler.readline()  # empty line

                            if target_sequence_end != end:
                                print(file)

                        # check sequence
                        allowed_bases = ["A", "U", "G", "C"]
                        new_structure = []
                        complete_seq = []
                        for base, base_str in zip(target_sequence, structure):
                            if base.upper() in allowed_bases:
                                new_structure.append(base_str)
                                complete_seq.append(base)

                        new_structure = "".join(new_structure)
                        complete_seq = "".join(complete_seq).upper()

                        add_loops, add_stem = find_loops(new_structure)
                        new_structure_shape = structure_description(add_loops, add_stem)

                        new_updated_structure = new_structure.replace('<', '('). \
                            replace('>', ')'). \
                            replace(',', '_'). \
                            replace('.', '_'). \
                            replace('-', '_')

                        if int(start) > int(end):
                            start, end = end, start

                        line_construct.extend((
                            candidate_id,
                            evalue, score,
                            start, end, strand,
                            structure, new_updated_structure,
                            new_structure_shape,
                            model_info,
                            target_sequence
                        ))

                        line_construct = '\t'.join(line_construct)
                        output.write(f"{line_construct}\n")

print("Done")
