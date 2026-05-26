"""Patch geny.py to support paired .gz FASTQ input."""

import sys


def patch(path: str) -> None:
    with open(path, "r") as f:
        content = f.read()

    # 1. Add --file2 argument to argparse
    old_argparse = (
        "    parser.add_argument('file', "
        'help="input BAM/CRAM/FASTA/FASTQ file", metavar=\'file\')'
    )
    new_argparse = (
        "    parser.add_argument('file', "
        'help="input BAM/CRAM/FASTA/FASTQ file", metavar=\'file\')\n'
        "    parser.add_argument('--file2', default=None, "
        'help="optional second FASTQ file (read 2) for paired-end input")'
    )
    assert old_argparse in content, "argparse target not found"
    content = content.replace(old_argparse, new_argparse)

    # 2. Replace get_fq_reads to support paired files
    old_fq = (
        "@timeit\n"
        "def get_fq_reads(path: str) -> Sample:\n"
        "    f = Sample(path)\n"
        "    log(f'[sample] fasta: {path}')\n"
        "    with pysam.FastxFile(path) as fq:\n"
        "        for read in fq:\n"
        '            name, pair = read.name.split("/")\n'
        "            f.reads.append(Sample.Read(\n"
        "                name,\n"
        '                pair == "2",\n'
        "                read.sequence,\n"
        "                read.comment,\n"
        "                '', # read.quality\n"
        "            ))\n"
        "    f.read_len = len(f.reads[0].seq)\n"
        "    f.reads.sort(key=lambda r: (r.name, r.pair))\n"
        "    return f"
    )
    new_fq = (
        "@timeit\n"
        "def get_fq_reads(path: str, path2: str = None) -> Sample:\n"
        "    f = Sample(path)\n"
        "    if path2:\n"
        "        log(f'[sample] paired fastq: {path} + {path2}')\n"
        "        for fq_path, is_read2 in [(path, False), (path2, True)]:\n"
        "            with pysam.FastxFile(fq_path) as fq:\n"
        "                for read in fq:\n"
        '                    name = read.name.split("/")[0]\n'
        "                    f.reads.append(Sample.Read(\n"
        "                        name,\n"
        "                        is_read2,\n"
        "                        read.sequence,\n"
        "                        read.comment,\n"
        "                        '',\n"
        "                    ))\n"
        "    else:\n"
        "        log(f'[sample] fasta: {path}')\n"
        "        with pysam.FastxFile(path) as fq:\n"
        "            for read in fq:\n"
        '                if "/" in read.name:\n'
        '                    name, pair = read.name.split("/")\n'
        '                    is_read2 = (pair == "2")\n'
        "                else:\n"
        "                    name = read.name\n"
        "                    is_read2 = False\n"
        "                f.reads.append(Sample.Read(\n"
        "                    name,\n"
        "                    is_read2,\n"
        "                    read.sequence,\n"
        "                    read.comment,\n"
        "                    '',\n"
        "                ))\n"
        "    f.read_len = len(f.reads[0].seq)\n"
        "    f.reads.sort(key=lambda r: (r.name, r.pair))\n"
        "    return f"
    )
    assert old_fq in content, "get_fq_reads target not found"
    content = content.replace(old_fq, new_fq)

    # 3. Fix extension check and pass file2
    old_main = (
        "    if path.endswith('.fa') or path.endswith('.fq') "
        "or path.endswith('.fastq'):\n"
        "        sample = get_fq_reads(path)\n"
        "        cov = args.coverage"
    )
    new_main = (
        "    fq_exts = ('.fa', '.fq', '.fastq', '.fa.gz', '.fq.gz', '.fastq.gz')\n"
        "    if path.endswith(fq_exts):\n"
        "        file2 = getattr(args, 'file2', None)\n"
        "        sample = get_fq_reads(path, file2)\n"
        "        cov = args.coverage"
    )
    assert old_main in content, "main extension check target not found"
    content = content.replace(old_main, new_main)

    with open(path, "w") as f:
        f.write(content)

    print("Patched successfully")


if __name__ == "__main__":
    patch(sys.argv[1])
