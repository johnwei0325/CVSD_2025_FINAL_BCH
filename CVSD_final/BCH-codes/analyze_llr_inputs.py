#!/usr/bin/env python3
#-*- coding: utf-8 -*-

"""
Script to analyze LLR test inputs from p100.txt and calculate error locations
for BCH code: m=6, n=63, k=51, t=2
"""

import sys
import argparse
from finitefield import *
from bch import *

# BCH code configurations
BCH_CONFIGS = {
    '100': {
        'name': '(63, 51) BCH',
        'm': 6,
        'n': 63,
        'k': 51,
        't': 2,
        'input_file': 'p100.txt',
        'answer_file': 'p100a.txt'
    },
    '200': {
        'name': '(255, 239) BCH',
        'm': 8,
        'n': 255,
        'k': 239,
        't': 2,
        'input_file': 'p200.txt',
        'answer_file': 'p200a.txt'
    },
    '300': {
        'name': '(1023, 983) BCH',
        'm': 10,
        'n': 1023,
        'k': 983,
        't': 4,
        'input_file': 'p300.txt',
        'answer_file': 'p300a.txt'
    }
}

def parse_8bit_signed_llr(llr_8bits):
    """
    Parse an 8-bit signed integer LLR value.
    
    :param llr_8bits: 8-bit binary string or integer
    :returns: Signed integer value (-128 to 127)
    """
    if isinstance(llr_8bits, str):
        value = int(llr_8bits, 2)
    else:
        value = llr_8bits
    
    # Convert to signed 8-bit integer (two's complement)
    if value & 0x80:  # MSB is set (negative)
        return value - 256
    else:
        return value

def llr_to_hard_decision(llr_value):
    """
    Convert LLR value to hard decision bit.
    - LLR > 0 (positive) -> favor bit 0
    - LLR < 0 (negative) -> favor bit 1
    - LLR == 0 -> favor bit 0 (default)
    
    :param llr_value: Signed integer LLR value
    :returns: Hard decision bit (0 or 1)
    """
    return 1 if llr_value < 0 else 0

def parse_llr_rows_to_codeword(llr_rows, n=63):
    """
    Parse LLR rows into a codeword.
    
    Each row contains 8 LLR values (8 bits each = 64 bits per row).
    LLR order: LLR0 (don't care), LLR1 (X^(n-1)), LLR2 (X^(n-2)), ..., LLR63 (X^0)
    
    :param llr_rows: List of rows, each row is 64 bits (8 LLR values)
    :param n: Codeword length
    :returns: Integer representing the hard decision codeword
    """
    # Concatenate all bits from rows
    all_bits = ''.join(llr_rows)
    
    # Parse into 8-bit LLR values
    llr_values = []
    for i in range(0, len(all_bits), 8):
        llr_8bits_str = all_bits[i:i+8]
        if len(llr_8bits_str) == 8:
            llr_signed = parse_8bit_signed_llr(llr_8bits_str)
            llr_values.append(llr_signed)
    
    # We need n+1 LLR values (LLR0 to LLRn)
    # LLR0 is don't care, so we use LLR1 to LLR63
    if len(llr_values) < n + 1:
        raise ValueError(f"Not enough LLR values: need {n+1}, got {len(llr_values)}")
    
    # Convert LLR values to hard decisions
    # LLR1 corresponds to X^(n-1) (MSB), LLR63 corresponds to X^0 (LSB)
    hard_decisions = []
    for i in range(1, n + 1):  # Skip LLR0, use LLR1 to LLR63
        hard_bit = llr_to_hard_decision(llr_values[i])
        hard_decisions.append(hard_bit)
    
    # Convert to integer: MSB first (LLR1) to LSB last (LLR63)
    codeword = 0
    for i, bit in enumerate(hard_decisions):
        codeword |= (bit << (n - 1 - i))  # LLR1 is MSB (position n-1), LLR63 is LSB (position 0)
    
    return codeword, llr_values

def get_error_locations(bch, received_codeword):
    """
    Calculate error locations for a received codeword without correcting it.
    
    :param bch: BCH code object
    :param received_codeword: Received codeword as integer
    :returns: List of error positions (0-indexed from LSB)
    """
    # Calculate syndromes
    syndromes, is_error = get_syndromes(
        primitive_polynomial=bch.primitive_polynomial,
        received_message=received_codeword,
        cyclotomic_cosets=bch.cyclotomic_cosets,
        logarithm_table=bch.logarithm_table,
        power=bch.power,
        t=bch.t,
        b=bch.b
    )
    
    if not is_error:
        return []  # No errors detected
    
    # Use Berlekamp-Massey to find error locator polynomial
    sigma = berlekamp_massey_decode(
        syndromes=syndromes,
        logarithm_table=bch.logarithm_table,
        power=bch.power,
        t=bch.t
    )
    
    # Find roots of error locator polynomial
    roots = find_roots_of_sigma(
        sigma=sigma,
        power=bch.power,
        logarithm_table=bch.logarithm_table
    )
    
    # Convert roots to error positions
    error_positions = get_error_positions(roots=roots, power=bch.power)
    
    return error_positions

def read_llr_file(filename):
    """
    Read LLR values from file.
    
    :param filename: Path to the LLR input file
    :returns: List of LLR strings (one per line)
    """
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f.readlines() if line.strip()]
    return lines

def read_answer_file(filename, t):
    """
    Read answer file containing expected error positions.
    
    Format: Each line is a binary representation of one error position.
    For each input, there are t lines (one per error position).
    
    :param filename: Path to the answer file
    :param t: Number of errors per input
    :returns: List of lists, where each inner list contains error positions for one input
    """
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f.readlines() if line.strip()]
    
    # Group lines by input (t lines per input)
    answers = []
    for i in range(0, len(lines), t):
        input_answers = []
        for j in range(i, min(i + t, len(lines))):
            answer_str = lines[j]
            try:
                # Each line is a binary number representing one error position
                error_pos = int(answer_str, 2)
                input_answers.append(error_pos)
            except ValueError:
                pass
        if input_answers:
            answers.append(input_answers)
    
    return answers

def extract_codewords_from_llr(llr_lines, n=63, num_inputs=2, mode='concatenate', verbose=False):
    """
    Extract codewords from LLR lines.
    
    Each row contains 8 LLR values (8-bit signed integers = 64 bits per row).
    For n=63, we need 64 LLR values (LLR0 to LLR63) = 8 rows.
    
    Mode 'concatenate': For two inputs with 16 rows total:
    - First 8 rows (lines 0-7) -> first input (64 LLR values)
    - Next 8 rows (lines 8-15) -> second input (64 LLR values)
    
    Mode 'per_row': Each row represents 8 LLR values (not used for n=63)
    
    :param llr_lines: List of LLR strings (64 bits each = 8 LLR values)
    :param n: Codeword length
    :param num_inputs: Number of inputs to extract (for concatenate mode)
    :param mode: 'concatenate' or 'per_row'
    :param verbose: If True, print detailed parsing information
    :returns: List of tuples (codeword, llr_values) where codeword is integer and llr_values is list
    """
    codewords_data = []
    
    if mode == 'per_row':
        # Each row has 8 LLR values, but we need 64 LLRs for n=63
        # This mode is not suitable for n=63, but kept for flexibility
        print("Warning: per_row mode may not work correctly for n=63")
        for line in llr_lines:
            if len(line) >= 64:
                codeword, llr_vals = parse_llr_rows_to_codeword([line], n=n)
                codewords_data.append((codeword, llr_vals))
    else:  # mode == 'concatenate'
        # Calculate rows per input based on n
        # Each row has 8 LLR values, we need (n+1) LLR values per input
        import math
        rows_per_input = math.ceil((n + 1) / 8)
        total_rows = len(llr_lines)
        
        for input_idx in range(num_inputs):
            start_row = input_idx * rows_per_input
            end_row = start_row + rows_per_input
            
            if verbose:
                print(f"\n{'='*70}")
                print(f"Input {input_idx + 1}: Reading rows {start_row} to {end_row-1} ({rows_per_input} rows with 64 bits each)")
                print(f"  Each row contains 8 LLR values (8-bit signed integers)")
                print(f"{'='*70}")
            
            # Show each row being read
            rows_for_input = llr_lines[start_row:end_row]
            if verbose:
                for i, row in enumerate(rows_for_input):
                    row_num = start_row + i
                    print(f"Row {row_num} (line {row_num+1}): {row} (length: {len(row)} bits)")
                    
                    # Show the 8 LLR values in this row
                    llr_in_row = []
                    for j in range(0, 64, 8):
                        llr_8bits = row[j:j+8]
                        llr_signed = parse_8bit_signed_llr(llr_8bits)
                        llr_in_row.append(llr_signed)
                    print(f"  LLR values in row: {llr_in_row}")
            
            # Parse LLR rows to codeword
            try:
                codeword, llr_values = parse_llr_rows_to_codeword(rows_for_input, n=n)
                codewords_data.append((codeword, llr_values))
                
                if verbose:
                    print(f"\nParsed {len(llr_values)} LLR values (LLR0 to LLR{len(llr_values)-1})")
                    print(f"  LLR0 (don't care): {llr_values[0]}")
                    print(f"  LLR1 to LLR{n} used for codeword")
                    print(f"\nCodeword (integer): {codeword}")
                    print(f"Codeword (binary, MSB to LSB): {bin(codeword)[2:].zfill(n)}")
                    
                    # Show hard decisions for first few and last few LLRs
                    print(f"\nHard decisions (LLR -> bit):")
                    for i in range(1, min(6, n+1)):
                        llr_val = llr_values[i]
                        hard_bit = llr_to_hard_decision(llr_val)
                        print(f"  LLR{i} = {llr_val:4d} (signed) -> bit {hard_bit} (X^{n-i})")
                    if n > 5:
                        print(f"  ...")
                        for i in range(max(1, n-4), n+1):
                            llr_val = llr_values[i]
                            hard_bit = llr_to_hard_decision(llr_val)
                            print(f"  LLR{i} = {llr_val:4d} (signed) -> bit {hard_bit} (X^{n-i})")
                        
            except ValueError as e:
                print(f"Error: {e}")
    
    return codewords_data

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Analyze LLR test inputs and calculate error locations for BCH codes',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
BCH Code Options:
  100  - (63, 51) BCH: m=6, n=63, k=51, t=2
  200  - (255, 239) BCH: m=8, n=255, k=239, t=2
  300  - (1023, 983) BCH: m=10, n=1023, k=983, t=4

Input/Output files:
  100  - p100.txt / p100a.txt
  200  - p200.txt / p200a.txt
  300  - p300.txt / p300a.txt
        """
    )
    parser.add_argument(
        'code',
        choices=['100', '200', '300'],
        help='BCH code selection: 100, 200, or 300'
    )
    parser.add_argument(
        '--testdata-dir',
        type=str,
        default='/Users/johnwei/Desktop/CVSD_final/1141_final/01_RTL/testdata',
        help='Directory containing test data files (default: 1141_final/01_RTL/testdata)'
    )
    parser.add_argument(
        '--compare',
        action='store_true',
        help='Compare results with answer file'
    )
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Show detailed LLR parsing information'
    )
    
    args = parser.parse_args()
    
    # Get BCH code configuration
    config = BCH_CONFIGS[args.code]
    m = config['m']
    n = config['n']
    k = config['k']
    t = config['t']
    dist = 2 * t + 1
    b = 1  # Starting power for cyclotomic cosets
    
    print("="*70)
    print(f"BCH Code: {config['name']}")
    print("="*70)
    print(f"Parameters: m = {m}, n = {n}, k = {k}, t = {t}")
    
    # Initialize BCH code
    try:
        primitive_poly = get_primitive_polynomial(m)
        bch = BCH(n, dist, b, primitive_poly)
        print(f"Primitive polynomial: {bin(primitive_poly)}")
    except Exception as e:
        print(f"Error initializing BCH code: {e}")
        return
    print()
    
    # Read LLR inputs
    input_file = f"{args.testdata_dir}/{config['input_file']}"
    print(f"Reading LLR inputs from: {input_file}")
    try:
        llr_lines = read_llr_file(input_file)
        print(f"Found {len(llr_lines)} lines of LLR data")
    except FileNotFoundError:
        print(f"Error: Input file not found: {input_file}")
        return
    except Exception as e:
        print(f"Error reading input file: {e}")
        return
    print()
    
    # Calculate rows per input
    # Each row has 64 bits = 8 LLR values (8 bits each)
    # We need (n+1) LLR values per input (LLR0 to LLRn)
    # So we need ceil((n+1)/8) rows per input
    import math
    rows_per_input = math.ceil((n + 1) / 8)
    total_rows_needed = rows_per_input * 2  # Assuming 2 inputs
    
    print(f"Rows per input: {rows_per_input} (need {n+1} LLR values, 8 per row)")
    print(f"Total rows needed for 2 inputs: {total_rows_needed}")
    print()
    
    # Extract codewords
    num_inputs = 2
    extraction_mode = 'concatenate'
    
    # Set verbose flag for extraction
    global_extraction_verbose = args.verbose
    
    codewords_data = extract_codewords_from_llr(llr_lines, n=n, num_inputs=num_inputs, mode=extraction_mode, verbose=args.verbose)
    print(f"\n{'='*70}")
    print(f"Summary: Extracted {len(codewords_data)} codeword(s) from {len(llr_lines)} rows (mode: {extraction_mode})")
    print(f"{'='*70}\n")
    
    # Read answer file if requested
    answer_data = None
    if args.compare:
        answer_file = f"{args.testdata_dir}/{config['answer_file']}"
        print(f"Reading expected answers from: {answer_file}")
        try:
            answer_data = read_answer_file(answer_file, t)
            print(f"Found {len(answer_data)} input(s) with error positions")
            for idx, ans in enumerate(answer_data):
                print(f"  Input {idx+1}: {sorted(ans)}")
        except FileNotFoundError:
            print(f"Warning: Answer file not found: {answer_file}")
        except Exception as e:
            print(f"Warning: Error reading answer file: {e}")
        print()
    
    # Process each codeword and find error locations
    for idx, (r, llr_vals) in enumerate(codewords_data):
        print(f"{'='*70}")
        print(f"Input {idx + 1}:")
        if args.verbose:
            print(f"  Received codeword r (binary): {bin(r)[2:].zfill(n)}")
            print(f"  Received codeword r (hex): {hex(r)}")
        print()
        
        # Calculate error locations
        error_positions = get_error_locations(bch, r)
        
        if len(error_positions) == 0:
            print("  No errors detected!")
        else:
            print(f"  Found {len(error_positions)} error(s) at position(s): {sorted(error_positions)}")
            print(f"  Error positions (0-indexed from LSB, MSB is position {n-1}):")
            for pos in sorted(error_positions):
                msb_pos = n - 1 - pos
                print(f"    - Bit position {pos} (LSB-indexed) / Position {msb_pos} (MSB-indexed)")
        
        # Compare with answer file if available
        if args.compare and answer_data and idx < len(answer_data):
            expected_positions = answer_data[idx]
            print(f"\n  Expected error positions (from answer file): {sorted(expected_positions)}")
            if set(error_positions) == set(expected_positions):
                print("  ✓ MATCH: Computed positions match expected positions!")
            else:
                print("  ✗ MISMATCH: Computed positions do not match expected positions")
                print(f"    Computed: {sorted(error_positions)}")
                print(f"    Expected: {sorted(expected_positions)}")
        
        # Also decode to see the corrected codeword
        if args.verbose:
            decoded_msg, corrected_codeword = bch.decode_ex(r)
            if r != corrected_codeword:
                print(f"\n  Corrected codeword: {bin(corrected_codeword)[2:].zfill(n)}")
                print(f"  Decoded message: {bin(decoded_msg)[2:].zfill(k)}")
        print()

if __name__ == '__main__':
    main()

