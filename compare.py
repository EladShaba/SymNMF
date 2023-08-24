import symnmf # your python symnmf script
import os
import numpy as np
import pandas as pd
import filecmp

def compare_directories(dir1, dir2):
  """Compares the two directories and returns a list of files that are different."""
  files1 = os.listdir(dir1)
  files2 = os.listdir(dir2)

  # Get the list of files that are in both directories
  common_files = list(set(files1) & set(files2))

  # Compare the files in the common_files list
  differences = []
  for file in common_files:
    file1 = os.path.join(dir1, file)
    file2 = os.path.join(dir2, file)
    if not filecmp.cmp(file1, file2, shallow=False):
      differences.append(file)

  return differences

def create_my_outputs():
    input_dir = "inputs"
    output_dir = "my_outputs"
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    
    # Open all the files in the input directory
    files = [os.path.join(input_dir, f) for f in os.listdir(input_dir)]

    # Run the symnmf algorithm on each input file
    for i, file in enumerate(files):
        max_iter = 300
        eps= 0.0001 # 1e-4
        df = pd.read_csv(file, header=None)
        n = df.shape[0]
        k = np.random.randint(1, 10)
        while k > n:
            k = np.random.randint(1, 10)
        data_points = df.to_numpy().tolist()
        # your function that returns the H_mat
        H_mat = symnmf.getSymnmf(data_points, n, k, eps, max_iter)
        
        # Write the output to a file
        with open(os.path.join(output_dir, f"outputs_{i}.txt"), "w") as f:
            for row in H_mat:
                line = ','.join([format(row[i], ".4f") for i in range(len(row))]) + "\n"
                f.write(line)
    
def main():
    np.random.seed(0)
    create_my_outputs()
    differences = compare_directories("outputs", "my_outputs")
    if differences != []:
        print(differences)
    else:
        print("all outputs are the same")


if __name__ == "__main__":
  main()