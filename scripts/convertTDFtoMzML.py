#!/usr/bin/env python

"""Conversion program to convert a Bruker TIMS file to a single mzML

"""
import argparse

from diapysef.convert_to_mzml import run


def main():

    parser = argparse.ArgumentParser(description ="Conversion program to convert a Bruker raw data file from a timsTOF Pro instrument into a single mzML.")
    parser.add_argument("-a", "--analysis_dir",
                        help = "The location of the directory containing raw data (usually .d)",
                        dest = 'analysis_dir',
                        required = True)
    parser.add_argument("-o", "--output_name",
                        help = "The name of the output file (mzML)",
                        dest = "output_fname",
                        required = True)
    parser.add_argument("-m", "--merge",
                        help = "Number of consecutive frames to sum up (squash). This is useful to boost S/N if exactly repeated frames are measured.",
                        type = int,
                        default = -1,
                        dest = "merge_scans")
    parser.add_argument("--keep_frames",
                        help = "Whether to store frames exactly as measured or split them into individual spectra by precursor isolation window (default is to split them - this is almost always what you want).",
                        type = bool,
                        default = False,
                        dest = "keep_frames")
    parser.add_argument("--verbose",
                        help = "Verbosity",
                        type = int,
                        default = -1,
                        dest = "verbosity")
    parser.add_argument("--overlap",
                        help = "How many overlapping windows were recorded for the same m/z window - will split the output into N output files.",
                        type = int,
                        default = -1,
                        dest = "overlap_scans")
    parser.add_argument("-r", "--framerange",
                        help = "The minimum and maximum Frames to convert. Useful to only convert a part of a file.",
                        type = int,
                        nargs = 2,
                        default = [-1, -1],
                        dest = "frame_limit")
    parser.add_argument("-z", "--use_zlib",
                        help = "Use zlib compression?",
                        type = bool,
                        default = False,
                        dest = "use_zlib")
    parser.add_argument("-n", "--use_numpress",
                        help = "Use numpress compression? NB - This uses hardcoded numpress params",
                        type = bool,
                        default = False,
                        dest = "use_numpress")
    args = parser.parse_args()
    run(**vars(args))


if __name__ == "__main__":
    main()

