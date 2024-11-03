# purpose of this script: annotate one mitogenome with set of CM models

import argparse
import subprocess
from pathlib import Path


def construct_cm(path_to_cm, path_to_output, path_to_fasta, cores):
    cmd = (
        f'cmsearch --max '  # turn off filtering
        f'--incE 0.05 '
        f'--notextw --smxsize 80000 '
        f'--cpu {cores} -g '  # global mode
        f'-o {path_to_output} '
        f'{path_to_cm} ',
        f'{path_to_fasta}'
    )
    return cmd


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, help='Path to mitogenome in fasta format', required=True)
    parser.add_argument('-o', '--output', type=str, help='Path to cmsearch output', required=True)
    parser.add_argument('-c', '--cms', type=str, help='Path to CMs set', required=True)
    parser.add_argument('-t', '--threads', type=str, help='Threads number for cmsearch', required=True)

    args = parser.parse_args()

    path_to_fasta = Path(args.input)
    path_to_output = Path(args.output)
    path_to_cms = Path(args.cms)
    threads_number = args.threads

    cms_batch = path_to_cms.glob("trn*.cm")

    for cm in cms_batch:

        cm_name = cm.name[:-3]
        print(cm_name)

        mitogenome_id = path_to_fasta.name.split(".fasta")[-1]

        # prepare output folder if needed
        if not path_to_output.exists():
            Path.mkdir(path_to_output)

        # check prev annotation
        output_annotation = Path(path_to_output, f"{cm_name}_annotation")

        if not output_annotation.exists():

            cmd_construct = construct_cm(cm,
                                         output_annotation,
                                         path_to_fasta,
                                         threads_number)
            cmd_construct = ' '.join(cmd_construct)

            try:
                process = subprocess.run(cmd_construct, check=True, shell=True)
            except subprocess.CalledProcessError as e:
                print(e)

    print("Done")
