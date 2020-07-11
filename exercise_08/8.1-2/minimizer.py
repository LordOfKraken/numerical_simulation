#!/usr/bin/env python

import os
def main():
    mu_min = 0.7
    mu_max = 1.0
    mu_step = 0.01
    mu_nstep = int((mu_max - mu_min)/mu_step)

    sigma_min = 0.5
    sigma_max = 0.8
    sigma_step = 0.01
    sigma_nstep = int((sigma_max - sigma_min)/sigma_step)

    x_start = 0

    for m in range(mu_nstep + 1):
        for s in range(sigma_nstep + 1):
            cmd = "./main.x {} {:.2f} {:.2f}".format(x_start,mu_min  + mu_step*m,sigma_min + sigma_step*s)
            print(cmd)
            os.system(cmd)

if __name__ == "__main__":
    main()
