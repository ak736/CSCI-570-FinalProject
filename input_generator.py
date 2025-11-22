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


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 2:
        print("Usage: python input_generator.py input.txt")
        sys.exit(1)

    X, Y = parse_input_file(sys.argv[1])

    print("\nGenerated String X:")
    print(X)
    print("Length X:", len(X))

    if Y:
        print("\nGenerated String Y:")
        print(Y)
        print("Length Y:", len(Y))
    else:
        print("\nWarning: Second string not found.")

