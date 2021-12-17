import os, sys
import matplotlib.pyplot as plt

if __name__ == '__main__':
    try:
        path = 'C:/Users/mokin/Desktop/git/CudaFEMproject/prak_results/'
        for filename in os.listdir(path):
            if filename.endswith('.txt'):
                with open(path + filename, 'r') as file_read:
                    data = file_read.readlines()
                data = [line.rstrip() for line in data]
                x, y = [], []
                for d in data:
                    x.append(float(d.split()[0]))
                    y.append(float(d.split()[1]))
                plt.plot(x, y)
        plt.show()
    except Exception as e:
        print('Something wrong with postprocessing')
        sys.exit(0)

    
