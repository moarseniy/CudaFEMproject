import os, sys
import matplotlib.pyplot as plt
import numpy as np

if __name__ == '__main__':
    leg_titles = []
    is_filter = True
    N = 60       # Number of entries for the moving average filter
    try:
        path = 'C:/Users/mexika/Documents/Qt_code/CudaFEMproject/prak_results/'
        for filename in os.listdir(path):
            if filename.endswith('.txt'):
                with open(path + filename, 'r') as file_read:
                    data = file_read.readlines()
                data = [line.rstrip() for line in data]
                x, y = [], []
                for d in data:
                    x.append(float(d.split()[0]))
                    y.append(float(d.split()[1]))
                if is_filter:
                    y_filtered = np.convolve(y, np.ones(N)/N, mode='valid')
                    if N % 2 == 0:
                        plt.plot(x[N//2-1:-N//2], y_filtered)
                    else:
                        plt.plot(x[(N-1)//2:-(N-1)//2], y_filtered)
                else:
                    plt.plot(x, y)
                leg_titles.append(filename)
        plt.legend(leg_titles)
        plt.show()
    except Exception as e:
        print('Something wrong with postprocessing')
        sys.exit(0)

    
