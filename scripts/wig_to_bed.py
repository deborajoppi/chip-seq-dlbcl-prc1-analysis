#!/usr/bin/env python3

import argparse
import gzip


def parse_args():
    parser = argparse.ArgumentParser(description="Convert variableStep WIG signal to merged BED intervals.")
    parser.add_argument("--input", required=True, help="Input .wig.gz file")
    parser.add_argument("--output", required=True, help="Output BED path")
    parser.add_argument("--min-signal", type=float, default=0.1, help="Minimum signal to retain")
    return parser.parse_args()


def main():
    args = parse_args()

    current_chrom = None
    current_span = None
    interval_chrom = None
    interval_start = None
    interval_end = None
    intervals_written = 0

    with gzip.open(args.input, "rt") as handle, open(args.output, "w") as out:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("track"):
                continue
            if line.startswith("variableStep"):
                fields = dict(part.split("=", 1) for part in line.split()[1:])
                current_chrom = fields["chrom"]
                current_span = int(fields.get("span", "1"))
                continue

            pos_str, value_str = line.split()
            pos = int(pos_str)
            value = float(value_str)

            if value < args.min_signal:
                continue

            start = pos - 1
            end = start + current_span

            if interval_chrom == current_chrom and interval_end == start:
                interval_end = end
            else:
                if interval_chrom is not None:
                    intervals_written += 1
                    out.write(f"{interval_chrom}\t{interval_start}\t{interval_end}\tH3K27me3_mask_{intervals_written}\t.\t.\n")
                interval_chrom = current_chrom
                interval_start = start
                interval_end = end

        if interval_chrom is not None:
            intervals_written += 1
            out.write(f"{interval_chrom}\t{interval_start}\t{interval_end}\tH3K27me3_mask_{intervals_written}\t.\t.\n")

    print(f"Wrote {intervals_written} merged H3K27me3 intervals to {args.output}")


if __name__ == "__main__":
    main()
