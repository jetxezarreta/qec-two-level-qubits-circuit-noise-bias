# Script aimed to understand bias loss in CNOT gate operation
import qutip as qu
import numpy as np
from scipy.io import savemat

# Define CNOT Hamiltonian from [Science 6,34 (2020)]

# ZX interaction Hamiltonian, standard implementation of a CNOT gate
# Hcnot = V [((I1+Z1)/2)\otimes I2 + ((I1-Z1)/2)\otimes X2]
V = 1 # normalized
Hcnot = V * (qu.tensor((qu.identity(2) + qu.sigmaz())/2,qu.identity(2)) + qu.tensor((qu.identity(2) - qu.sigmaz())/2,qu.sigmax()))

# the CNOT gate occurs at VT = pi/2
T = 0.5 * np.pi / V

# Now we define the dissipators of the master equation
# We will study a biased model rateZ >> rate X
# by lambdax we refer to the non pure Z errors

eta = 100
lambdaz = eta * 0.002 / (3*(1+eta))
lambdax = 0.002 / (12*(1+eta))
#bias = lambdaz / lambdax

# define collapse operators
# Biased elements
Lopsz1 = qu.tensor(np.sqrt(lambdaz) * qu.sigmaz(), qu.identity(2))
Lopsz2 = qu.tensor(qu.identity(2), np.sqrt(lambdaz) * qu.sigmaz())
Lopsz3 = qu.tensor(qu.sigmaz(), np.sqrt(lambdaz) * qu.sigmaz())

# non-biased elements
Lopsx1 = qu.tensor(np.sqrt(lambdax) * qu.sigmax(), qu.identity(2))
Lopsx2 = qu.tensor(qu.identity(2), np.sqrt(lambdax) * qu.sigmax())
Lopsx3 = qu.tensor(qu.sigmax(), np.sqrt(lambdax) * qu.sigmax())
Lopsx4 = qu.tensor(qu.sigmay(), np.sqrt(lambdax) * qu.sigmax())
Lopsx5 = qu.tensor(qu.sigmaz(), np.sqrt(lambdax) * qu.sigmax())
Lopsx6 = qu.tensor(qu.identity(2), np.sqrt(lambdax) * qu.sigmay())
Lopsx7 = qu.tensor(np.sqrt(lambdax) * qu.sigmay(), qu.identity(2))
Lopsx8 = qu.tensor(np.sqrt(lambdax) * qu.sigmax(), qu.sigmay())
Lopsx9 = qu.tensor(np.sqrt(lambdax) * qu.sigmay(), qu.sigmay())
Lopsx10 = qu.tensor(np.sqrt(lambdax) * qu.sigmax(), qu.sigmaz())
Lopsx11 = qu.tensor(np.sqrt(lambdax) * qu.sigmay(), qu.sigmaz())
Lopsx12 = qu.tensor(np.sqrt(lambdax) * qu.sigmaz(), qu.sigmay())


# Pauli transfer matrix
# it inputs a rho equal to a pauli matrix and then the obsevable is the othe Pauli matrix
# make this in a loop to evaluate the whole transfer matrix
outPaul = ([qu.tensor(qu.identity(2),qu.identity(2)), qu.tensor(qu.identity(2),qu.sigmax()),
             qu.tensor(qu.identity(2),qu.sigmay()), qu.tensor(qu.identity(2),qu.sigmaz()),
             qu.tensor(qu.sigmax(),qu.identity(2)), qu.tensor(qu.sigmax(),qu.sigmax()),
             qu.tensor(qu.sigmax(),qu.sigmay()), qu.tensor(qu.sigmax(),qu.sigmaz()),
             qu.tensor(qu.sigmay(),qu.identity(2)), qu.tensor(qu.sigmay(),qu.sigmax()),
             qu.tensor(qu.sigmay(),qu.sigmay()), qu.tensor(qu.sigmay(),qu.sigmaz()),
             qu.tensor(qu.sigmaz(),qu.identity(2)), qu.tensor(qu.sigmaz(),qu.sigmax()),
             qu.tensor(qu.sigmaz(),qu.sigmay()),qu.tensor(qu.sigmaz(),qu.sigmaz())])
inPaul = ([qu.tensor(qu.identity(2),qu.identity(2)), qu.tensor(qu.identity(2),qu.sigmax()),
             qu.tensor(qu.identity(2),qu.sigmay()), qu.tensor(qu.identity(2),qu.sigmaz()),
             qu.tensor(qu.sigmax(),qu.identity(2)), qu.tensor(qu.sigmax(),qu.sigmax()),
             qu.tensor(qu.sigmax(),qu.sigmay()), qu.tensor(qu.sigmax(),qu.sigmaz()),
             qu.tensor(qu.sigmay(),qu.identity(2)), qu.tensor(qu.sigmay(),qu.sigmax()),
             qu.tensor(qu.sigmay(),qu.sigmay()), qu.tensor(qu.sigmay(),qu.sigmaz()),
             qu.tensor(qu.sigmaz(),qu.identity(2)), qu.tensor(qu.sigmaz(),qu.sigmax()),
             qu.tensor(qu.sigmaz(),qu.sigmay()),qu.tensor(qu.sigmaz(),qu.sigmaz())])
transfMatNoise = np.zeros((16,16))
idx = 0
jdx = 0
for inp in inPaul:
    jdx = 0
    for outp in outPaul:
        result = qu.mesolve(
                    Hcnot, inp, [0, T], [Lopsz1, Lopsz2, Lopsz3, Lopsx1, Lopsx2, Lopsx3, Lopsx4, Lopsx5, Lopsx6, Lopsx7, Lopsx8, Lopsx9, Lopsx10, Lopsx11, Lopsx12], outp
                )
        Rij = result.expect
        transfMatNoise[jdx,idx] = np.multiply(Rij[0][1],1/4)
        jdx += 1
    idx += 1
print(transfMatNoise)



# Result ideal (Pauli transfer matrix)
transfMatIdeal = np.zeros((16,16))
idx = 0
jdx = 0
for inp in inPaul:
    jdx = 0
    for outp in outPaul:
        result = qu.mesolve(
                    Hcnot, inp, [0, T], [], [outp]
                )
        Rij = result.expect
        transfMatIdeal[jdx,idx] = np.multiply(Rij[0][1],1/4)
        jdx += 1
    idx += 1
print(transfMatNoise)

# save data
savemat("PauliTransfMatCNOTtestIdeal.mat", {'IdealRij': transfMatIdeal})
savemat("PauliTransfMatCNOTtestNoisy.mat", {'NoisyRij': transfMatNoise})

