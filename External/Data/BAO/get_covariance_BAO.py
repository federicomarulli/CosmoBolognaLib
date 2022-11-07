import numpy as np

data = np.genfromtxt("BAO_Addison2013.dat", unpack=True, usecols=(0, 1, 2))

ndata = len(data[0])
corr, cov = np.zeros(shape=(ndata, ndata)), np.zeros(shape=(ndata, ndata))


for i in range(ndata):
    corr[i][i]=1

# SDSS DR7

corr[1][2] = 0.57
corr[2][1] = corr[1][2]

#WiggleZ

WiggleZ_inv_cov = [ [24532.1, -25137.7, 12099.1],
                    [-25137.7, 134598.4, -64783.9],
                    [12099.1, -64783.9, 128837.6]]

WiggleZ_cov = np.linalg.inv(WiggleZ_inv_cov)

for i in range(len(WiggleZ_cov)):
    print np.sqrt(WiggleZ_cov[i][i]), data[2][i+3]
    for j in range(len(WiggleZ_cov[i])):
        corr[i+3][j+3] = WiggleZ_cov[i][j]/np.sqrt(WiggleZ_cov[i][i]*WiggleZ_cov[j][j])

#BOSS DR9

corr[6][7] = 0.55
corr[7][6] = corr[6][7]

print corr

ff = open("BAO_Addison2013_covariance.dat", "w")
ff.write("#redshift   redshift   covariance   correlation \n")

#Full correlation
for i in range(ndata):
    for j in range(ndata):
        cov[i][j] = data[2][i]*data[2][j]*corr[i][j]
        ll = "%2.10f  %2.10f  %2.10f  %2.10f \n"%(data[0][i], data[0][j], data[2][i]*data[2][j]*corr[i][j], corr[i][j])
        ff.write(ll)
    ff.write("\n")

ff.close()
print np.linalg.inv(cov)   

