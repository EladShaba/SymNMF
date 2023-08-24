import analysis
import os
import numpy as np
import pandas as pd


def main():
    np.random.seed(0)
    
    input_dir = "inputs"
    output_dir = "outputs"
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    
    # Open all the files in the input directory
    files = [os.path.join(input_dir, f) for f in os.listdir(input_dir)]

    # Run the symnmf algorithm on each file
    for i, file in enumerate(files):
        max_iter = 300
        eps= 0.0001 # 1e-4
        df = pd.read_csv(file, header=None)
        n = df.shape[0]
        k = np.random.randint(1, 10)
        while k > n:
            k = np.random.randint(1, 10)
        data_points = df.to_numpy().tolist()
        H_mat = analysis.getSymnmf(data_points, n, k, eps, max_iter)

        # Write the output to a file
        with open(os.path.join(output_dir, f"outputs_{i}.txt"), "w") as f:
            for row in H_mat:
                line = ','.join([format(row[i], ".4f") for i in range(len(row))]) + "\n"
                f.write(line)

if __name__ == "__main__":
  main()