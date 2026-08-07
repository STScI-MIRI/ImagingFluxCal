import matplotlib.pyplot as plt
import numpy as np


if __name__ == "__main__":

    # make plot
    fontsize = 16
    font = {"size": fontsize}
    plt.rc("font", **font)
    plt.rc("lines", linewidth=2)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("ytick.major", width=2)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))

    loss = 0.03
    x = np.arange(20.0)
    y = np.power(1.0 - loss, x)

    ax.plot(x, y, "k-")
    ax.plot(x, (1.0 - loss * x), "b--")

    plt.show()