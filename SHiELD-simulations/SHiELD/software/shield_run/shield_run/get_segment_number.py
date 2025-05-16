import argparse

from shield_run.append import get_segment_number
from pathlib import Path


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("destination", type=Path)
    args, extra_args = parser.parse_known_args()

    segment_number = get_segment_number(args.destination)
    print(segment_number)
