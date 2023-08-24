import os
import numpy as np

def write_matrices_to_files(matrices, output_dir):
    """Writes the given matrices to .txt files in the given output directory."""
    for i, matrix in enumerate(matrices):
        with open(os.path.join(output_dir, f"input_{i}.txt"), "w") as f:
            for row in matrix:
                line = ','.join([format(row[i], "f") for i in range(len(row))]) + "\n"
                f.write(line)

def main():
    # Generate some matrices
    shapes = [(np.random.randint(2, 500), np.random.randint(2, 20)) for i in range(50)]
    matrices = [np.random.uniform(-10, 10, (r, c)) for r, c in shapes]
        
    # Create the output directory
    output_dir = "inputs"
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # Write the matrices to files
    write_matrices_to_files(matrices, output_dir)

if __name__ == "__main__":
    np.random.seed(0)
    main()