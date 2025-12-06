import sys
import time
import psutil


# Mismatch costs alpha_pq
# Rows/Cols: A, C, G, T. Order is fixed to map to indices 0, 1, 2, 3.
CHAR_TO_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

#     A   C   G   T
ALPHA = [
    [0, 110, 48, 94],  # A
    [110, 0, 118, 48],  # C
    [48, 118, 0, 110],  # G
    [94, 48, 110, 0]   # T
]

# Gap penalty delta
DELTA = 30


def generate_string(base, indices):
    s = base
    for idx in indices:
        s = s[:idx + 1] + s + s[idx + 1:]
    return s


def parse_input_file(file_path):
    with open(file_path, "r") as f:
        lines = [line.strip() for line in f if line.strip() != ""]

    i = 0

    # First base string
    s0 = lines[i]
    i += 1

    # Read indices for first string
    s_indices = []
    while i < len(lines) and lines[i].isdigit():
        s_indices.append(int(lines[i]))
        i += 1

    X = generate_string(s0, s_indices)

    # Second base string
    if i >= len(lines):
        return X, None

    t0 = lines[i]
    i += 1

    # Read indices for second string
    t_indices = []
    while i < len(lines) and lines[i].isdigit():
        t_indices.append(int(lines[i]))
        i += 1

    Y = generate_string(t0, t_indices)

    return X, Y


def get_mismatch_cost(char1, char2):
    """Looks up the mismatch cost alpha_pq for two characters."""
    if char1 == '-' or char2 == '-':
        return DELTA

    idx1 = CHAR_TO_INDEX.get(char1)
    idx2 = CHAR_TO_INDEX.get(char2)

    if idx1 is not None and idx2 is not None:
        return ALPHA[idx1][idx2]

    # Should not happen based on problem constraints
    return float('inf')


# --- BASIC DP ALGORITHM (Needleman-Wunsch-like) ---

def basic_sequence_alignment(S, T):
    """
    Calculates the minimum alignment cost and reconstructs the optimal alignment
    using standard Dynamic Programming (O(m*n) time and space).
    """

    m = len(S)
    n = len(T)

    # OPT[i][j] = min cost of aligning S[0..i-1] with T[0..j-1]
    OPT = [[0] * (n + 1) for _ in range(m + 1)]

    # Initialize column 0 and row 0
    for i in range(1, m + 1):
        OPT[i][0] = OPT[i-1][0] + DELTA

    for j in range(1, n + 1):
        OPT[0][j] = OPT[0][j-1] + DELTA

    # Fill DP Table
    for i in range(1, m + 1):
        for j in range(1, n + 1):

            # 1. Mismatch/Match: OPT(i-1, j-1) + alpha(x_i, y_j)
            cost_match_mismatch = OPT[i-1][j-1] + \
                get_mismatch_cost(S[i-1], T[j-1])

            # 2. Gap in T (x_i aligns with gap): OPT(i-1, j) + delta
            cost_gap_T = OPT[i-1][j] + DELTA

            # 3. Gap in S (y_j aligns with gap): OPT(i, j-1) + delta
            cost_gap_S = OPT[i][j-1] + DELTA

            OPT[i][j] = min(cost_match_mismatch, cost_gap_T, cost_gap_S)

    min_cost = OPT[m][n]

    # Traceback to reconstruct the optimal alignment
    aligned_S = []
    aligned_T = []
    i, j = m, n

    while i > 0 or j > 0:
        current_cost = OPT[i][j]

        if i > 0 and j > 0:
            # Check for match/mismatch move (Diagonal)
            match_mismatch_cost = OPT[i-1][j-1] + \
                get_mismatch_cost(S[i-1], T[j-1])

            if current_cost == match_mismatch_cost:
                aligned_S.append(S[i-1])
                aligned_T.append(T[j-1])
                i -= 1
                j -= 1
                continue

        if i > 0:
            # Check for gap in T
            gap_T_cost = OPT[i-1][j] + DELTA

            if current_cost == gap_T_cost:
                aligned_S.append(S[i-1])
                aligned_T.append('_')
                i -= 1
                continue

        if j > 0:
            # Check for gap in S
            gap_S_cost = OPT[i][j-1] + DELTA

            if current_cost == gap_S_cost:
                aligned_S.append('_')
                aligned_T.append(T[j-1])
                j -= 1
                continue

    # The traceback builds the alignment in reverse order, so reverse the results
    return min_cost, "".join(aligned_S[::-1]), "".join(aligned_T[::-1])


# --- TIME AND MEMORY TRACKING / MAIN EXECUTION ---

def process_memory():
    """Get current memory usage in kilobytes."""
    process = psutil.Process()
    memory_info = process.memory_info()
    return int(memory_info.rss / 1024)


def main(input_filepath, output_filepath):

    # Measure memory before algorithm
    mem_before = process_memory()

    # Time
    start_time = time.time() * 1000  # Time in milliseconds

    try:
        # Generate Input Strings (S and T)
        string_S, string_T = parse_input_file(input_filepath)

        # Run Basic DP Algorithm
        cost, aligned_S, aligned_T = basic_sequence_alignment(
            string_S, string_T)

        end_time = time.time() * 1000
        time_taken = end_time - start_time

        # Measure memory after algorithm
        mem_after = process_memory()

        # Calculate memory used by algorithm
        # memory_consumed_kb = mem_after - mem_before
        memory_consumed_kb = mem_after  # Total memory used at end

    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        return

    # --- Write Output File ---
    try:
        with open(output_filepath, 'w') as outfile:
            # 1. Cost of the alignment (Integer)
            outfile.write(f"{cost}\n")
            # 2. First string alignment
            outfile.write(f"{aligned_S}\n")
            # 3. Second string alignment
            outfile.write(f"{aligned_T}\n")
            # 4. Time in Milliseconds (Float)
            outfile.write(f"{time_taken:.3f}\n")
            # 5. Memory in Kilobytes (Float)
            outfile.write(f"{memory_consumed_kb}\n")
    except Exception as e:
        print(f"Error writing to output file: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(
            "Usage: python3 basic.py <input_file_path> <output_file_path>", file=sys.stderr)
        sys.exit(1)

    input_file_path = sys.argv[1]
    output_file_path = sys.argv[2]

    main(input_file_path, output_file_path)
