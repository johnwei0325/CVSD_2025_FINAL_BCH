# BCH LLR Input Analyzer

This script analyzes LLR (Log-Likelihood Ratio) test inputs and calculates error locations for BCH codes.

## Usage

```bash
python3 analyze_llr_inputs.py <code> [options]
```

## BCH Code Options

- `100` - (63, 51) BCH: m=6, n=63, k=51, t=2
  - Input: `p100.txt`
  - Answer: `p100a.txt`

- `200` - (255, 239) BCH: m=8, n=255, k=239, t=2
  - Input: `p200.txt`
  - Answer: `p200a.txt`

- `300` - (1023, 983) BCH: m=10, n=1023, k=983, t=4
  - Input: `p300.txt`
  - Answer: `p300a.txt`

## Options

- `--compare` - Compare results with answer file
- `--verbose` - Show detailed LLR parsing information
- `--testdata-dir <path>` - Specify test data directory (default: `1141_final/01_RTL/testdata`)

## Examples

```bash
# Basic usage for (63,51) BCH
python3 analyze_llr_inputs.py 100

# With answer comparison
python3 analyze_llr_inputs.py 100 --compare

# With verbose output
python3 analyze_llr_inputs.py 100 --verbose

# For (255,239) BCH with comparison
python3 analyze_llr_inputs.py 200 --compare

# For (1023,983) BCH with comparison
python3 analyze_llr_inputs.py 300 --compare
```

## Output

The script will:
1. Parse LLR values from the input file (8-bit signed integers)
2. Convert LLR values to hard decisions (received codeword)
3. Calculate error locations using BCH decoding
4. Optionally compare with expected answers from the answer file

## Requirements

- Python 3
- Files from `BCH-codes/` directory: `bch.py`, `finitefield.py`

