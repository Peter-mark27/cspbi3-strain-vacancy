"""Validate the included data and redraw selected manuscript figures."""
import argparse
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def figure_numbers(value):
    selected = set()
    try:
        for part in value.split(','):
            if '-' in part:
                start, end = map(int, part.split('-', 1))
                if start > end:
                    raise ValueError
                selected.update(range(start, end + 1))
            else:
                selected.add(int(part))
        if not selected or not selected.issubset(range(1, 11)):
            raise ValueError
    except ValueError as exc:
        raise argparse.ArgumentTypeError('Use figure numbers 1-10, for example 3-10 or 1,3,5.') from exc
    return sorted(selected)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--figures', type=figure_numbers, default=figure_numbers('1-10'))
    parser.add_argument('--output-dir', type=Path, default=ROOT / 'reproduced_figures')
    args = parser.parse_args()
    output = args.output_dir.resolve()
    if output == (ROOT / 'figures').resolve():
        parser.error('Choose a different output directory to preserve the manuscript reference PDFs.')
    output.mkdir(parents=True, exist_ok=True)
    environment = os.environ.copy()
    environment['MPLCONFIGDIR'] = str(output / '.matplotlib')
    environment['PYTHONDONTWRITEBYTECODE'] = '1'
    subprocess.run([sys.executable, '-B', str(ROOT / 'scripts' / 'validate_data.py')], check=True, env=environment)
    for number in args.figures:
        print(f'Redrawing figure {number}...', flush=True)
        subprocess.run([sys.executable, '-B', str(ROOT / 'scripts' / f'figure{number}.py'),
                        '--output-dir', str(output)], check=True, env=environment)
    print(f'Completed {len(args.figures)} figures. Results: {output}', flush=True)


if __name__ == '__main__':
    main()
