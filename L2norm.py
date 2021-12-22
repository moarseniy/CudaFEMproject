import os, sys
import numpy as np

def Simpson(func, dp):
    N = len(func)
    integral_sum = 0.0
    
    if N % 2 == 1:
        k_end = (N-1)//2
    else:
        k_end = N//2 - 1
    
    for k in range(k_end):
        integral_sum += func[2*k] + 4.0*func[2*k + 1] + func[2*k + 2]                        
    
    integral_sum = dp/3.0*integral_sum
    
    if N % 2 == 0:
        integral_sum += dp/2.0*(func[N - 2]+func[N - 1])

    return integral_sum

def L2norm(x, w):
    if len(x) != len(w):
        print('Mesh x and function don''t match!')
        exit()
    
    func = np.power(np.abs(w), 2)
    
    h = x[1] - x[0]                # spacing
           
    return np.sqrt(Simpson(func, h))


if __name__ == '__main__':
    is_filter = False
    N = 150       # Number of entries for the moving average filter
    funcs = {}
    Idx = 1
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
                    y = np.convolve(y, np.ones(N)/N, mode='valid')
                    if N % 2 == 0:
                        x = x[N//2-1:-N//2]
                    else:
                        x = x[(N-1)//2:-(N-1)//2]

                funcs[Idx] = (x,y)
                print(Idx,'\t',filename)
                Idx += 1
    except Exception as e:
        print('Something wrong with postprocessing')
        sys.exit(0)

    for i in range(1, Idx - 1):
        x = funcs[i][0]
        y1 = np.array(funcs[i][1])
        y2 = np.array(funcs[i+1][1])
        norm_y2 = L2norm(x,y2)
        norm_diff = L2norm(x,y1-y2)
        print('L2 norm (y',i,', y',i+1,') =\t', L2norm(x,y1-y2),' ( ', norm_diff/norm_y2*100,'% )\t L2 norm (y',i+1,') =\t',norm_y2, sep='')

    
