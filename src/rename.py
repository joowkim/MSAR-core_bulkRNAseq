import argparse
import logging
import os
import re
from collections import OrderedDict, defaultdict
from dataclasses import dataclass
from typing import List, Dict, DefaultDict

logging.basicConfig(level=logging.DEBUG,
                    format=' %(asctime)s - %(levelname)s- %(message)s')


def read_sample_sheet(sample_sheet: str) -> OrderedDict[str, str]:
    """

    :param sample_sheet:
    :return: result_dict: key - sample_name, value = sample_id
    """
    result_dict: OrderedDict[str, str] = OrderedDict()
    assert os.path.isfile(sample_sheet), f"{sample_sheet} is not found!"

    with open(sample_sheet) as fin:
        flag: bool = False  # this is for not including the columns - Sample_ID,Sample_Name ...
        for line in fin:
            if "Sample_ID" in line:
                flag: bool = True
                # don't include the columns - Sample_ID,Sample_Name...
                line: str = fin.readline()

            if flag:
                # line looks like this
                # Sample_ID,Sample_Name,Sample_Plate,Sample_Well,Index_Plate_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
                # WGS1,GC-LX-1434_1,,,,,TAAGGCGA,,AGAGGATA,GC-LX-1434,
                tmp_line_list: List[str] = line.split(",")
                sample_id: str = tmp_line_list[0]
                sample_name: str = tmp_line_list[1]

                if sample_id in result_dict:
                    print(f"duplicated sample_id found - {sample_id} in {sample_sheet}")
                    raise ValueError
                elif sample_name in result_dict:
                    print(f"duplicated sample_name found - {sample_name} in {sample_sheet}")
                    raise ValueError
                else:
                    result_dict[sample_name] = sample_id
    return result_dict


def get_fq_path(fastq_dir: str) -> Dict[str, List[str]]:
    """

    :param fastq_dir:
    :return: result_dict: key - file_name, value - relative_file_path
    """
    result_dict: DefaultDict[str] = defaultdict(list)

    assert os.path.isdir(fastq_dir), f"{fastq_dir} is not found!"

    for parent_dir, sub_dir, files in os.walk(fastq_dir):
        for fq in files:
            if fq.endswith("fastq.gz"):
                abs_path: str = os.path.join(parent_dir, fq)
                assert os.path.isfile(abs_path), f"{abs_path} is not valid!"

                # get the prefix of the name of fastq files
                file_name: str = re.sub("_S[\d]+_L00[1-4]_R[1-3]_00[1-2].fastq.gz", "", fq)

                if file_name.endswith("fastq.gz"):
                    file_name: str = re.sub("_R[1-3].fastq.gz", "", fq)
                    result_dict[file_name].append(abs_path)
                else:
                    result_dict[file_name].append(abs_path)

    return result_dict


def generate_rename_fq_dict(sample_sheet_dict: Dict, fq_path_dict: Dict) -> OrderedDict[str]:
    """

    :param sample_sheet_dict: OrderedDict([('GC-LX-1434_1', 'WGS1'), ('GC-LX-1434_2', 'WGS2'), sample_name is key; sample_id is value.
    :param fq_path_dict: defaultdict(<class 'list'>, {'GC.LX.1434_31': ['/Users/kimj32/tmp/change_fq/GC.LX.1434/Sample_GC.LX.1434_31/GC.LX.1434_31_S31_L002_R2_001.fastq.gz', '/Users/kimj32/tmp/change_fq/GC.LX.1434/Sample_GC.LX.1434_31/GC.LX.1434_31_S31_L001_R1_001.fastq.gz'...
    :return:
    """
    result_dict: OrderedDict = OrderedDict()

    for sample_name in sample_sheet_dict:
        repl_sample_name: str = sample_name.replace("-", ".")

        if not repl_sample_name in fq_path_dict:
            print(f"{sample_name} is not found in the name of fastq files.")
            # raise ValueError

        if repl_sample_name in fq_path_dict:
            # new_sample_name: str = f"{sample_sheet_dict.get(sample_name, 'NA')}-{repl_sample_name}"
            new_sample_name: str = f"{sample_sheet_dict.get(sample_name, 'NA')}"
            # logging.info(f"new_sample_name - {new_sample_name}")
            fq_file_list: List[str] = fq_path_dict.get(repl_sample_name)

            for fq_path in fq_file_list:
                fq_dir_name: str = os.path.dirname(fq_path)
                assert os.path.isdir(fq_dir_name), f"{fq_dir_name} is not found!"
                fq_basename: str = os.path.basename(fq_path)
                new_fq_basename: str = fq_basename.replace(repl_sample_name, new_sample_name)
                # logging.info(f"fq_basename - {fq_basename}")
                # logging.info(f"new_fq_basename - {new_fq_basename}")
                final_new_fq_path: str = os.path.join(fq_dir_name, new_fq_basename)
                # logging.info(f"ori_fq_path - {fq_path}")
                # logging.info(f"new_fq_path - {final_new_fq_path}")

                result_dict[fq_path] = final_new_fq_path

    return result_dict


def rename(fq_path_dict: Dict, dry_run: bool, fq_path_for_log_file: str):
    """

    :param fq_path_dict: key - original name of fq, value - new name of fq
    :param dry_run:
    :return:
    """

    log_file_prefix: str = "rename.log"
    log_file_path:str = os.path.join(fq_path_for_log_file, log_file_prefix)

    if dry_run:
        for ori_path in fq_path_dict:
            new_path: str = fq_path_dict.get(ori_path)
            logging.info(f"this is a dry run. rename {ori_path} -> {new_path}.")
    else:
        with open(log_file_path, 'w') as fout:
            fout.write('original_name,new_name\n')
            for ori_path in fq_path_dict:
                new_path: str = fq_path_dict.get(ori_path)
                logging.info(f"rename {ori_path} -> {new_path}")
                os.rename(ori_path, new_path)
                fout.write(f"{ori_path},{new_path}\n")


def roll_back(log_file: str):
    assert os.path.isfile(log_file), f"f{log_file} is not found."

    with open(log_file)as fin:
        fin.readline()
        for line in fin:
            line: str = line.strip()
            tmp_list: List[str] = line.split(",")
            ori_path: str = tmp_list[0]
            new_path: str = tmp_list[1]
            os.rename(new_path, ori_path)
    print("roll-back is done.")


@dataclass
class Args:
    sample_sheet: str
    fq_path: str
    dry_run: bool
    roll_back: bool
    log_file: str


def get_args() -> Args:
    parser = argparse.ArgumentParser(
        description="rename fastq files - [sample_id]_[sample_name]...fastq.gz",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        "-s",
        "--sample_sheet",
        help="sample_sheet.csv",
        metavar="sample_sheet",
        type=str,
        required=True,
    )

    parser.add_argument(
        "-fq",
        "--fq_path",
        help="fastq_path",
        metavar="fq_path",
        type=str,
        required=True,
    )

    parser.add_argument(
        "-d",
        "--dry_run",
        help="dry_run",
        action="store_true",
        default=False,
    )

    parser.add_argument(
        "-rb",
        "--roll_back",
        help="roll back",
        action="store_true",
        default=False,
    )

    parser.add_argument(
        '-l',
        "--log_file",
        help='log_file for roll back',
        type=str,
        default="",
    )

    args = parser.parse_args()
    return Args(args.sample_sheet, args.fq_path, args.dry_run, args.roll_back, args.log_file)


def main():
    args = get_args()
    sample_sheet: str = args.sample_sheet
    fq_path: str = args.fq_path
    dry_run: bool = args.dry_run
    roll_back_flag: bool = args.roll_back
    log_file: str = args.log_file

    if roll_back_flag:
        roll_back(log_file=log_file)
    else:
        sample_sheet_dict: Dict = read_sample_sheet(sample_sheet=sample_sheet)
        fq_dict: Dict = get_fq_path(fastq_dir=fq_path)
        rename_fq_dict: Dict = generate_rename_fq_dict(sample_sheet_dict=sample_sheet_dict, fq_path_dict=fq_dict)
        rename(fq_path_dict=rename_fq_dict, dry_run=dry_run, fq_path_for_log_file=fq_path)


if __name__ == '__main__':
    main()

## useful command for checking file names
## cat rename_log.log| sed 's@/Users/kimj32/tmp/change_fq/GC.LX.1434/Sample_GC.LX.1434_[0-9]*/@@g' | sed 's@[a-zA-Z0-9_]*-GC.LX@GC.LX@g' > sed.csv