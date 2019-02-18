import numpy
from matplotlib import pyplot as plt

if __name__ == '__main__':
    x = [1, 2, 3]
    plt.ion() # turn on interactive mode
    for loop in range(0,3):
        y = numpy.dot(x, loop)
        plt.figure()
        plt.plot(x,y)
        plt.show()
        _ = raw_input("Press [enter] to continue.")
