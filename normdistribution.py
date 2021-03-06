

cent = [29.657,27.924,-0.15,20.242,17.552,18.35,7.499,7.062,25.895,15.21,
15.776,13.489,11.912,8.417,3.761,10.252,9.339,4.351,13.223,7.49,7.361,11.351,
14.798,11.695,6.382,11.262,5.369,6.437,11.884,5.634,2.967,11.774,5.897,2.237,
2.695,2.318,2.95,3.134,1.009,-1.837,1.941,-3.279,7.18,-0.068,-0.17,1.14,-2.796,
2.791,5.428,9.542,-8.951,-14.331
]

cl15=[-3.216,3.666,0.916,-9.27,-11.233,19.664,4.715,-2.951,21.801,0.672,3.239,-12.307,
6.789,4.799,-8.826,-5.909,-3.735,-21.165,-14.167,0.969,4.728,9.597,6.377,-11.435,-3.182,
-3.011,4.538,-5.628,-8.004,8.088,7.802,-13.896,-16.128,-7.757,9.683,-17.505,-18.874,-2.574,
-1.62,-17.5,-20.7,-7.613,8.944,4.304,-10.647,-8.506,-8.133,-13.447,0.844,-15.969,-13.032,6.947
]
cl85=[8.33,7.63,-1.005,-2.082,-8.522,16.666,-3.195,6.424,32.39,2.849,4.331,-11.733,4.203,5.642,
-15.395,-2.179,-6.578,-22.035,-8.425,0.226,9.401,10.169,14.765,-3.199,-2.411,-8.897,0.509,-2.27,
-4.709,-4.535,8.602,-3.421,-4.782,-13.572,4.41,-17.258,-21.517,2.174,-2.863,-18.674,-11.717,-9.456,
3.082,-4.056,-8.632,-8.706,-5.247,-14.918,3.331,-12.221,-9.574,1.758
]

import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt


if __name__ == "__main__": 
    ## cent, cl15, cl85
    x = sorted(cent) 

    x = sorted(cl15+cl85) 

    # Generate some data for this demonstration.
    data = x

    # Fit a normal distribution to the data:
    mu, std = norm.fit(data)

    # Plot the histogram.
    plt.hist(data, bins=25, density=True, alpha=0.6, color='gray')

    # Plot the PDF.
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    p = norm.pdf(x, mu, std)
    plt.plot(x, p, 'k', linewidth=2)
    title = "Fit results: mu = %.2f,  std = %.2f" % (mu, std)
    plt.title(title)

    plt.show()