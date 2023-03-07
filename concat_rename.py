import argparse
import logging
import os
import re
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from typing import List, Dict, DefaultDict

logging.basicConfig(
    format="%(asctime)s,%(msecs)d %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s",
    datefmt="%Y-%m-%d:%H:%M:%S",
    level=logging.DEBUG,
)


def get_fq_path(fastq_dir: str) -> Dict[str, List[str]]:
    """
    :param fastq_dir:
    :return: result_dict: key - file_name, value - relative_file_path
    """
    result_dict: DefaultDict[str, List[str]] = defaultdict(list)

    assert os.path.isdir(fastq_dir), f"{fastq_dir} is not found!"

    for parent_dir, sub_dir, files in os.walk(fastq_dir):
        for fq in files:
            if fq.endswith("fastq.gz"):
                abs_path: str = os.path.join(parent_dir, fq)
                assert os.path.isfile(abs_path), f"{abs_path} is not valid!"

                # get the prefix of the name of fastq files
                file_name: str = re.sub("_S[\d]+_L00[1-2]_R[1-2]_00[1-2].fastq.gz", "", fq)
                result_dict[file_name].append(abs_path)

    return result_dict


@dataclass
class FastqInfo:
    sample_name: str
    submit_name: str
    r1_list: List[str]
    r2_list: List[str]


def get_r1_r2_fastq(fastq_path_dict: Dict, sample_key_dict: Dict) -> List[FastqInfo]:
    result_list: List[FastqInfo] = list()
    # sample_name from iLab
    for sample_name in sorted(fastq_path_dict):
        if sample_name in sample_key_dict:
            fq_list: List[str] = fastq_path_dict.get(sample_name)
            r1_list: List[str] = sorted([i for i in fq_list if "L001_R1" in i or "L002_R1" in i])
            r2_list: List[str] = sorted([i for i in fq_list if "L001_R2" in i or "L002_R2" in i])
            ins: FastqInfo = FastqInfo(sample_name=sample_name,
                                       submit_name=sample_key_dict.get(sample_name),
                                       r1_list=r1_list,
                                       r2_list=r2_list,
                                       )
            result_list.append(ins)
        else:
            print(f"{sample_name} is not found in the sample_key.csv")
            raise ValueError

    return result_list


def concat_fastq(fastq_r1_r2_list: List[FastqInfo]):
    os.makedirs("rename", exist_ok=True)
    os.makedirs("logs/rename", exist_ok=True)

    for fastq_info_ins in fastq_r1_r2_list:
        r1_list: List[str] = fastq_info_ins.r1_list
        r1_l1_fastq: str = r1_list[0]
        r1_l2_fastq: str = r1_list[1]
        r1_final: str = os.path.join("rename", f"{fastq_info_ins.submit_name}_R1.fastq.gz")
        # r2_list: List[str] = fastq_info_ins.r2_list
        cmd: str = f"cat {r1_l1_fastq} {r1_l2_fastq} > {r1_final}"
        logging.info(cmd)
        sub_out = subprocess.run([cmd], capture_output=True, shell=True)

        r2_list: List[str] = fastq_info_ins.r2_list
        r2_l1_fastq: str = r2_list[0]
        r2_l2_fastq: str = r2_list[1]
        r2_final: str = os.path.join("rename", f"{fastq_info_ins.submit_name}_R2.fastq.gz")
        cmd: str = f"cat {r2_l1_fastq} {r2_l2_fastq} > {r2_final}"
        logging.info(cmd)
        sub_out = subprocess.run([cmd], capture_output=True, shell=True)


def read_sample_key(csv_file: str) -> Dict[str, str]:
    assert os.path.isfile(csv_file), f"{csv_file} is not found."

    result_dict: Dict[str, str] = dict()

    with open(csv_file) as fin:
        line_flag: bool = False
        for line in fin:
            if line.startswith("Submitted_Name"):
                line_flag: bool = True
            if line_flag:
                # line looks like this
                # header
                # Submitted_Name, Sample_Name, index, index2
                # Nic_3,GC.NG.1429_6,TATGTAGTCA,CATTAGTGCG
                line: str = line.strip()
                tmp_line_list: List[str] = line.split(",")
                submit_name: str = tmp_line_list[0]
                # iLab Name
                sample_name: str = tmp_line_list[1]
                result_dict[sample_name] = submit_name
    return result_dict


def test():
    csv_file = "GC.NG.1429_SampleKey.csv"
    sample_key_dict = read_sample_key(csv_file)
    fq_path_dict = get_fq_path("rawdata")
    fq_list = get_r1_r2_fastq(fastq_path_dict=fq_path_dict, sample_key_dict=sample_key_dict)
    concat_fastq(fq_list)


@dataclass
class Args:
    sample_sheet: str
    raw_data: str


def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description="merge fastq if same libraries split into different lanes - Lane1/2",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-s",
                        "--sample_sheet",
                        help="sample_key.csv",
                        metavar="file",
                        type=str,
                        required=True
                        )

    parser.add_argument("-r",
                        "--raw_data",
                        help="raw_data_dir",
                        metavar="raw_data",
                        type=str,
                        required=True
                        )
    args = parser.parse_args()

    return Args(args.sample_sheet, args.raw_data)


def main():
    args = get_args()
    sample_key: str = args.sample_sheet
    raw_data: str = args.raw_data
    sample_key_dict = read_sample_key(sample_key)
    fq_path_dict = get_fq_path(raw_data)
    fq_list = get_r1_r2_fastq(fastq_path_dict=fq_path_dict, sample_key_dict=sample_key_dict)
    concat_fastq(fq_list)


if __name__ == '__main__':
    main()
    # test()
