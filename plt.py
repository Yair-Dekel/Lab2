import matplotlib.pyplot as plt

def read_numbers_from_file(file_path):
    with open(file_path, 'r') as file:
        numbers = [float(line.strip()) for line in file]
    return numbers

def plot_numbers(numbers, file_path):
    plt.figure(figsize=(10, 5))
    plt.plot(numbers, marker='o', linestyle='-', color='b')
    #plt.title('distance over 80 generations SUDOKU')
    plt.title(file_path + ' 200 generations')
    plt.xlabel('Generation')
    plt.ylabel('Value')
    plt.grid(True)
    plt.show()

# SUDOKU
# BINS
#this code takes a file of data and plots it 

if __name__ == "__main__":
    file_path = 'average_gene_distance_sudoku.txt'  # path to your file
    numbers = read_numbers_from_file(file_path)
    plot_numbers(numbers, file_path)
    file_path = 'variance_BINS.txt'  # path to your file
    numbers = read_numbers_from_file(file_path)
    plot_numbers(numbers, file_path)
    file_path = 'average_BINS.txt'  # path to your file
    numbers = read_numbers_from_file(file_path)
    plot_numbers(numbers, file_path)
    file_path = 'allels_BINS.txt'  # path to your file
    numbers = read_numbers_from_file(file_path)
    plot_numbers(numbers, file_path)
    file_path = 'fitness_BINS.txt'  # path to your file
    numbers = read_numbers_from_file(file_path)
    plot_numbers(numbers, file_path)    

    



