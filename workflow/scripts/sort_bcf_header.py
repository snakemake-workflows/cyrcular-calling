with open(snakemake.input.header, "rt") as header_file:
    header = [l.strip() for l in header_file.readlines()]
    file_format_line = header[0]
    chrom_line = header[-1]
    other_lines = header[1:-1]
    kinds = [
        "fileDate",
        "source",
        "reference",
        "contig",
        "phasing",
        "FILTER",
        "INFO",
        "FORMAT",
        "ALT",
        "assembly",
        "META",
        "SAMPLE",
        "PEDIGREE",
        "pedigreeDB",
    ]
    categories = {kind: [] for kind in kinds}
    others = []
    for line in other_lines:
        if "=" in line:
            kind = line.split("=")[0].lstrip("#")
            group = categories.get(kind, others)
        else:
            group = others
        group.append(line)

    with open(snakemake.output.sorted_header, "wt") as out:
        print(file_format_line, file=out)
        for kind in kinds:
            lines = categories[kind]
            for line in lines:
                print(line, file=out)
        for line in others:
            print(line, file=out)
        print(chrom_line, file=out)
