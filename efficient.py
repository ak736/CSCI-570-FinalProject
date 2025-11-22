#!/usr/bin/env python3
"""
Sequence Alignment using Hirschberg's Algorithm
Memory-efficient implementation using Divide and Conquer
CSCI-570 Fall 2025 Final Project
"""

import sys
import time
import psutil
import os

# ============================================================================
# CONSTANTS - HARDCODED AS PER PROJECT REQUIREMENTS
# ============================================================================

GAP_PENALTY = 30

# Mismatch cost matrix for DNA nucleotides (A, C, G, T)
MISMATCH_COSTS = {
    ('A', 'A'): 0,   ('A', 'C'): 110, ('A', 'G'): 48,  ('A', 'T'): 94,
    ('C', 'A'): 110, ('C', 'C'): 0,   ('C', 'G'): 118, ('C', 'T'): 48,
    ('G', 'A'): 48,  ('G', 'C'): 118, ('G', 'G'): 0,   ('G', 'T'): 110,
    ('T', 'A'): 94,  ('T', 'C'): 48,  ('T', 'G'): 110, ('T', 'T'): 0
}

# ============================================================================
# INPUT GENERATION (from specification)
# ============================================================================


def generate_string(base, indices):
    """
    Generate string by iteratively inserting the string into itself
    at specified indices (0-indexed).

    Args:
        base: Base string to start with
        indices: List of indices where to insert

    Returns:
        Final generated string after all insertions
    """
    current = base
    for idx in indices:
        # Insert current string into itself after position idx
        current = current[:idx + 1] + current + current[idx + 1:]
    return current


def parse_input_file(filepath):
    """
    Parse input file to generate two DNA strings for alignment.

    File format:
        base_string_1
        index_1
        index_2
        ...
        base_string_2
        index_1
        index_2
        ...

    Returns:
        Tuple of (string1, string2)
    """
    with open(filepath, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    idx = 0

    # Parse first string
    base1 = lines[idx]
    idx += 1

    indices1 = []
    while idx < len(lines) and lines[idx].replace('-', '').isdigit():
        indices1.append(int(lines[idx]))
        idx += 1

    string1 = generate_string(base1, indices1)

    # Parse second string
    if idx >= len(lines):
        return string1, ""

    base2 = lines[idx]
    idx += 1

    indices2 = []
    while idx < len(lines) and lines[idx].replace('-', '').isdigit():
        indices2.append(int(lines[idx]))
        idx += 1

    string2 = generate_string(base2, indices2)

    return string1, string2

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================


def get_mismatch_cost(char1, char2):
    """Get the mismatch cost for aligning two characters."""
    return MISMATCH_COSTS.get((char1, char2), 0)


def process_memory():
    """Get current memory usage in kilobytes."""
    process = psutil.Process()
    memory_info = process.memory_info()
    return int(memory_info.rss / 1024)

# ============================================================================
# SPACE-EFFICIENT DP FUNCTIONS (Core of Hirschberg's Algorithm)
# ============================================================================


def compute_alignment_score_forward(X, Y):
    """
    Compute the last row of the DP table using O(n) space.
    This represents the cost of aligning X with all prefixes of Y.

    Args:
        X: First string (length m)
        Y: Second string (length n)

    Returns:
        Array of length n+1 with alignment costs
    """
    m, n = len(X), len(Y)

    # We only need two rows: previous and current
    prev_row = [j * GAP_PENALTY for j in range(n + 1)]
    curr_row = [0] * (n + 1)

    for i in range(1, m + 1):
        # First column: cost of aligning X[0:i] with empty string
        curr_row[0] = i * GAP_PENALTY

        for j in range(1, n + 1):
            # Three options:
            # 1. Match/mismatch X[i-1] with Y[j-1]
            match_cost = prev_row[j - 1] + \
                get_mismatch_cost(X[i - 1], Y[j - 1])
            # 2. Gap in Y (delete from X)
            delete_cost = prev_row[j] + GAP_PENALTY
            # 3. Gap in X (insert into X)
            insert_cost = curr_row[j - 1] + GAP_PENALTY

            curr_row[j] = min(match_cost, delete_cost, insert_cost)

        # Swap rows for next iteration
        prev_row, curr_row = curr_row, prev_row

    return prev_row


def compute_alignment_score_backward(X, Y):
    """
    Compute alignment scores for suffixes by reversing both strings
    and running forward DP.

    Args:
        X: First string
        Y: Second string

    Returns:
        Array with costs for aligning suffixes, properly reversed
    """
    # Reverse both strings
    X_rev = X[::-1]
    Y_rev = Y[::-1]

    # Run forward DP on reversed strings
    scores = compute_alignment_score_forward(X_rev, Y_rev)

    # Reverse the result to get correct suffix costs
    return scores[::-1]

# ============================================================================
# HIRSCHBERG'S DIVIDE AND CONQUER ALGORITHM
# ============================================================================


def find_optimal_split(X, Y):
    """
    Find the optimal position to split Y when dividing the problem.

    Args:
        X: First string (will be split in middle)
        Y: Second string (split position to be found)

    Returns:
        Index in Y where optimal split occurs
    """
    m = len(X)
    mid = m // 2

    # Split X at midpoint
    X_left = X[:mid]
    X_right = X[mid:]

    # Compute forward scores: aligning X_left with prefixes of Y
    forward_scores = compute_alignment_score_forward(X_left, Y)

    # Compute backward scores: aligning X_right with suffixes of Y
    backward_scores = compute_alignment_score_backward(X_right, Y)

    # Find the split in Y that minimizes total cost
    n = len(Y)
    min_cost = float('inf')
    split_position = 0

    for j in range(n + 1):
        total_cost = forward_scores[j] + backward_scores[j]
        if total_cost < min_cost:
            min_cost = total_cost
            split_position = j

    return split_position


def hirschberg_recursive(X, Y):
    """
    Hirschberg's algorithm: recursively compute optimal alignment
    using divide and conquer with O(m+n) space.

    Args:
        X: First DNA string
        Y: Second DNA string

    Returns:
        Tuple of (aligned_X, aligned_Y)
    """
    m, n = len(X), len(Y)

    # BASE CASE 1: Empty X
    if m == 0:
        return '_' * n, Y

    # BASE CASE 2: Empty Y
    if n == 0:
        return X, '_' * m

    # BASE CASE 3: Single character in X
    if m == 1:
        # Try aligning X[0] with each position in Y
        min_cost = float('inf')
        best_alignment = None

        # Option: Align X[0] with some Y[j], gaps elsewhere
        for j in range(n):
            cost = (j * GAP_PENALTY +
                    get_mismatch_cost(X[0], Y[j]) +
                    (n - j - 1) * GAP_PENALTY)
            if cost < min_cost:
                min_cost = cost
                aligned_X = '_' * j + X[0] + '_' * (n - j - 1)
                aligned_Y = Y
                best_alignment = (aligned_X, aligned_Y)

        # Option: All gaps (X[0] at start, rest gaps)
        cost = GAP_PENALTY + n * GAP_PENALTY
        if cost < min_cost:
            best_alignment = (X[0] + '_' * n, '_' + Y)

        return best_alignment

    # BASE CASE 4: Single character in Y (optimization)
    if n == 1:
        min_cost = float('inf')
        best_alignment = None

        for i in range(m):
            cost = (i * GAP_PENALTY +
                    get_mismatch_cost(X[i], Y[0]) +
                    (m - i - 1) * GAP_PENALTY)
            if cost < min_cost:
                min_cost = cost
                aligned_X = X
                aligned_Y = '_' * i + Y[0] + '_' * (m - i - 1)
                best_alignment = (aligned_X, aligned_Y)

        # Option: All gaps (Y[0] at start, rest gaps)
        cost = m * GAP_PENALTY + GAP_PENALTY
        if cost < min_cost:
            best_alignment = ('_' + X, Y[0] + '_' * m)

        return best_alignment

    # RECURSIVE CASE: Divide and conquer
    mid = m // 2
    split_y = find_optimal_split(X, Y)

    # Divide both strings at optimal positions
    X_left = X[:mid]
    X_right = X[mid:]
    Y_left = Y[:split_y]
    Y_right = Y[split_y:]

    # Recursively solve left subproblem
    left_X, left_Y = hirschberg_recursive(X_left, Y_left)

    # Recursively solve right subproblem
    right_X, right_Y = hirschberg_recursive(X_right, Y_right)

    # Combine results
    aligned_X = left_X + right_X
    aligned_Y = left_Y + right_Y

    return aligned_X, aligned_Y


def calculate_alignment_cost(aligned_X, aligned_Y):
    """
    Calculate the total cost of an alignment.
    Used to get the final cost value.

    Args:
        aligned_X: First aligned string (with gaps '_')
        aligned_Y: Second aligned string (with gaps '_')

    Returns:
        Total alignment cost
    """
    if len(aligned_X) != len(aligned_Y):
        raise ValueError("Aligned strings must have equal length")

    total_cost = 0
    for i in range(len(aligned_X)):
        char_x = aligned_X[i]
        char_y = aligned_Y[i]

        if char_x == '_' and char_y == '_':
            raise ValueError("Both positions cannot be gaps")
        elif char_x == '_':
            total_cost += GAP_PENALTY
        elif char_y == '_':
            total_cost += GAP_PENALTY
        else:
            total_cost += get_mismatch_cost(char_x, char_y)

    return total_cost

# ============================================================================
# OUTPUT FUNCTIONS
# ============================================================================


def write_output(filepath, cost, aligned_X, aligned_Y, time_ms, memory_kb):
    """
    Write results to output file in required format.
    Creates directories if they don't exist.

    Format:
        Line 1: Alignment cost (integer)
        Line 2: First aligned string
        Line 3: Second aligned string
        Line 4: Time in milliseconds (float)
        Line 5: Memory in kilobytes (float)
    """
    # Create directory if it doesn't exist
    output_dir = os.path.dirname(filepath)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with open(filepath, 'w') as f:
        f.write(f"{int(cost)}\n")
        f.write(f"{aligned_X}\n")
        f.write(f"{aligned_Y}\n")
        f.write(f"{time_ms}\n")
        f.write(f"{memory_kb}\n")

# ============================================================================
# MAIN EXECUTION
# ============================================================================


def main():
    """Main function to run the memory-efficient sequence alignment."""

    # Check command-line arguments
    if len(sys.argv) != 3:
        sys.stderr.write(
            "Usage: python3 efficient.py <input_file> <output_file>\n")
        sys.exit(1)

    input_filepath = sys.argv[1]
    output_filepath = sys.argv[2]

    # Parse input and generate strings
    string_X, string_Y = parse_input_file(input_filepath)

    # Measure time and run algorithm
    start_time = time.time()
    aligned_X, aligned_Y = hirschberg_recursive(string_X, string_Y)
    end_time = time.time()

    # Measure memory after algorithm completes
    memory_used = process_memory()

    # Calculate metrics
    time_taken = (end_time - start_time) * 1000  # Convert to milliseconds

    # Calculate alignment cost
    alignment_cost = calculate_alignment_cost(aligned_X, aligned_Y)

    # Write output to file
    write_output(output_filepath, alignment_cost, aligned_X, aligned_Y,
                 time_taken, memory_used)


if __name__ == "__main__":
    main()
