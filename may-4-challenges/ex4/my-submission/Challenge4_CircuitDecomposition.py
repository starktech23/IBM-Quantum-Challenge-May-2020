#!/usr/bin/env python
# coding: utf-8

# # Exercise 4: Circuit Decomposition
# Wow! If you managed to solve the first three exercises, congratulations! The fourth problem is supposed to puzzle even the quantum experts among you, so don’t worry if you cannot solve it. If you can, hats off to you!
# 
# You may recall from your quantum mechanics course that quantum theory is unitary. Therefore, the evolution of any (closed) system can be described by a unitary. But given an arbitrary unitary, can you actually implement it on your quantum computer?
# 
# **"A set of quantum gates is said to be universal if any unitary transformation of the quantum data can be efficiently approximated arbitrarily well as a sequence of gates in the set."** (https://qiskit.org/textbook/ch-algorithms/defining-quantum-circuits.html)
# 
# Every gate you run on the IBM Quantum Experience is transpiled into single qubit rotations and CNOT (CX) gates. We know that these constitute a universal gate set, which implies that any unitary can be implemented using only these gates. However, in general it is not easy to find a good decomposition for an arbitrary unitary. Your task is to find such a decomposition.
# 
# You are given the following unitary:

# In[95]:


from may4_challenge.ex4 import get_unitary
U = get_unitary()
print(U)
print("U has shape", U.shape)


# #### What circuit would make such a complicated unitary?
# 
# Is there some symmetry, or is it random? We just updated Qiskit with the introduction of a quantum circuit library (https://github.com/Qiskit/qiskit-terra/tree/master/qiskit/circuit/library). This library gives users access to a rich set of well-studied circuit families, instances of which can be used as benchmarks (quantum volume), as building blocks in building more complex circuits (adders), or as tools to explore quantum computational advantage over classical computation (instantaneous quantum polynomial complexity circuits).

# In[2]:


from qiskit import QuantumCircuit
from may4_challenge.ex4 import check_circuit, submit_circuit


# **Using only single qubit rotations and CNOT gates, find a quantum circuit that approximates that unitary $U$ by a unitary $V$ up to an error $\varepsilon = 0.01$, such that $\lVert U - V\rVert_2 \leq \varepsilon$ !** 
# 
# Note that the norm we are using here is the spectral norm, $\qquad \lVert A \rVert_2 = \max_{\lVert \psi \rVert_2= 1} \lVert A \psi \rVert$.
# 
# This can be seen as the largest scaling factor that the matrix $A$ has on any initial (normalized) state $\psi$. One can show that this norm corresponds to the largest singular value of $A$, i.e., the square root of the largest eigenvalue of the matrix $A^\dagger A$, where $A^{\dagger}$ denotes the conjugate transpose of $A$.
# 
# **When you submit a circuit, we remove the global phase of the corresponding unitary $V$ before comparing it with $U$ using the spectral norm. For example, if you submit a circuit that generates $V = \text{e}^{i\theta}U$, we remove the global phase $\text{e}^{i\theta}$ from $V$ before computing the norm, and you will have a successful submission. As a result, you do not have to worry about matching the desired unitary, $U$, up to a global phase.**
# 
# As the single-qubit gates have a much higher fidelity than the two-qubit gates, we will look at the number of CNOT-gates, $n_{cx}$, and the number of u3-gates, $n_{u3}$, to determine the cost of your decomposition as 
# 
# $$
# \qquad \text{cost} = 10 \cdot n_{cx} + n_{u3}
# $$
# 
# Try to optimize the cost of your decomposition. 
# 
# **Note that you will need to ensure that your circuit is composed only of $u3$ and $cx$ gates. The exercise is considered correctly solved if your cost is smaller than 1600.**
# 
# ---
# For useful tips to complete this exercise as well as pointers for communicating with other participants and asking questions, please take a look at the following [repository](https://github.com/qiskit-community/may4_challenge_exercises). You will also find a copy of these exercises, so feel free to edit and experiment with these notebooks.
# 
# ---

# In[164]:


##### build your quantum circuit here
qc = QuantumCircuit(4)

#import qiskit.quantum_info as qi
#unitary = qi.Operator(qc)
#print(unitary)

unitary_H = [[ 0.25+0.j,  0.25+0.j,  0.25+0.j,  0.25+0.j,  0.25+0.j,
            0.25+0.j,  0.25+0.j,  0.25+0.j,  0.25+0.j,  0.25+0.j,
            0.25+0.j,  0.25+0.j,  0.25+0.j,  0.25+0.j,  0.25+0.j,
            0.25+0.j],
          [ 0.25+0.j, -0.25+0.j,  0.25+0.j, -0.25+0.j,  0.25+0.j,
           -0.25+0.j,  0.25+0.j, -0.25+0.j,  0.25+0.j, -0.25+0.j,
            0.25+0.j, -0.25+0.j,  0.25+0.j, -0.25+0.j,  0.25+0.j,
           -0.25+0.j],
          [ 0.25+0.j,  0.25+0.j, -0.25+0.j, -0.25+0.j,  0.25+0.j,
            0.25+0.j, -0.25+0.j, -0.25+0.j,  0.25+0.j,  0.25+0.j,
           -0.25+0.j, -0.25+0.j,  0.25+0.j,  0.25+0.j, -0.25+0.j,
           -0.25+0.j],
          [ 0.25+0.j, -0.25+0.j, -0.25+0.j,  0.25+0.j,  0.25+0.j,
           -0.25+0.j, -0.25+0.j,  0.25+0.j,  0.25+0.j, -0.25+0.j,
           -0.25+0.j,  0.25+0.j,  0.25+0.j, -0.25+0.j, -0.25+0.j,
            0.25+0.j],
          [ 0.25+0.j,  0.25+0.j,  0.25+0.j,  0.25+0.j, -0.25+0.j,
           -0.25+0.j, -0.25+0.j, -0.25+0.j,  0.25+0.j,  0.25+0.j,
            0.25+0.j,  0.25+0.j, -0.25+0.j, -0.25+0.j, -0.25+0.j,
           -0.25+0.j],
          [ 0.25+0.j, -0.25+0.j,  0.25+0.j, -0.25+0.j, -0.25+0.j,
            0.25+0.j, -0.25+0.j,  0.25+0.j,  0.25+0.j, -0.25+0.j,
            0.25+0.j, -0.25+0.j, -0.25+0.j,  0.25+0.j, -0.25+0.j,
            0.25+0.j],
          [ 0.25+0.j,  0.25+0.j, -0.25+0.j, -0.25+0.j, -0.25+0.j,
           -0.25+0.j,  0.25+0.j,  0.25+0.j,  0.25+0.j,  0.25+0.j,
           -0.25+0.j, -0.25+0.j, -0.25+0.j, -0.25+0.j,  0.25+0.j,
            0.25+0.j],
          [ 0.25+0.j, -0.25+0.j, -0.25+0.j,  0.25+0.j, -0.25+0.j,
            0.25+0.j,  0.25+0.j, -0.25+0.j,  0.25+0.j, -0.25+0.j,
           -0.25+0.j,  0.25+0.j, -0.25+0.j,  0.25+0.j,  0.25+0.j,
           -0.25+0.j],
          [ 0.25+0.j,  0.25+0.j,  0.25+0.j,  0.25+0.j,  0.25+0.j,
            0.25+0.j,  0.25+0.j,  0.25+0.j, -0.25+0.j, -0.25+0.j,
           -0.25+0.j, -0.25+0.j, -0.25+0.j, -0.25+0.j, -0.25+0.j,
           -0.25+0.j],
          [ 0.25+0.j, -0.25+0.j,  0.25+0.j, -0.25+0.j,  0.25+0.j,
           -0.25+0.j,  0.25+0.j, -0.25+0.j, -0.25+0.j,  0.25+0.j,
           -0.25+0.j,  0.25+0.j, -0.25+0.j,  0.25+0.j, -0.25+0.j,
            0.25+0.j],
          [ 0.25+0.j,  0.25+0.j, -0.25+0.j, -0.25+0.j,  0.25+0.j,
            0.25+0.j, -0.25+0.j, -0.25+0.j, -0.25+0.j, -0.25+0.j,
            0.25+0.j,  0.25+0.j, -0.25+0.j, -0.25+0.j,  0.25+0.j,
            0.25+0.j],
          [ 0.25+0.j, -0.25+0.j, -0.25+0.j,  0.25+0.j,  0.25+0.j,
           -0.25+0.j, -0.25+0.j,  0.25+0.j, -0.25+0.j,  0.25+0.j,
            0.25+0.j, -0.25+0.j, -0.25+0.j,  0.25+0.j,  0.25+0.j,
           -0.25+0.j],
          [ 0.25+0.j,  0.25+0.j,  0.25+0.j,  0.25+0.j, -0.25+0.j,
           -0.25+0.j, -0.25+0.j, -0.25+0.j, -0.25+0.j, -0.25+0.j,
           -0.25+0.j, -0.25+0.j,  0.25+0.j,  0.25+0.j,  0.25+0.j,
            0.25+0.j],
          [ 0.25+0.j, -0.25+0.j,  0.25+0.j, -0.25+0.j, -0.25+0.j,
            0.25+0.j, -0.25+0.j,  0.25+0.j, -0.25+0.j,  0.25+0.j,
           -0.25+0.j,  0.25+0.j,  0.25+0.j, -0.25+0.j,  0.25+0.j,
           -0.25+0.j],
          [ 0.25+0.j,  0.25+0.j, -0.25+0.j, -0.25+0.j, -0.25+0.j,
           -0.25+0.j,  0.25+0.j,  0.25+0.j, -0.25+0.j, -0.25+0.j,
            0.25+0.j,  0.25+0.j,  0.25+0.j,  0.25+0.j, -0.25+0.j,
           -0.25+0.j],
          [ 0.25+0.j, -0.25+0.j, -0.25+0.j,  0.25+0.j, -0.25+0.j,
            0.25+0.j,  0.25+0.j, -0.25+0.j, -0.25+0.j,  0.25+0.j,
            0.25+0.j, -0.25+0.j,  0.25+0.j, -0.25+0.j, -0.25+0.j,
            0.25+0.j]]


import numpy as np
S = np.transpose(unitary_H)
#I = np.matmul(S, unitary_H)
#print(I)

#V = np.matmul(T, S)
#T = np.matmul(unitary_H, U)


#V1 = V.tolist()
qc.h(0)
qc.h(1)
qc.h(2)
qc.h(3)

U1 = U.diagonal()
U1 = U1.tolist()
qc.diagonal(U1, [0, 1, 2, 3])

qc.h(0)
qc.h(1)
qc.h(2)
qc.h(3)

#V = U.tolist()
#qc.diagonal(V, [0, 1, 2, 3])
#V = [-0.21338834764831813+0.33838834764831804j, -0.21338834764831813+0.33838834764831804j, -0.21338834764831813+0.33838834764831804j, -0.21338834764831813+0.33838834764831804j, -0.21338834764831813+0.33838834764831804j, -0.21338834764831813+0.33838834764831804j, -0.21338834764831813+0.33838834764831804j, -0.21338834764831813+0.33838834764831804j, -0.21338834764831813+0.33838834764831804j, -0.21338834764831813+0.33838834764831804j, -0.21338834764831813+0.33838834764831804j, -0.21338834764831813+0.33838834764831804j, -0.21338834764831813+0.33838834764831804j, -0.21338834764831813+0.33838834764831804j, -0.21338834764831813+0.33838834764831804j, -0.21338834764831813+0.33838834764831804j]
#qc.diagonal(V, [0, 1, 2, 3])

# apply operations to your quantum circuit here

from qiskit.compiler import transpile
qc = transpile(qc, basis_gates = ['u3', 'cx'], optimization_level=3)
# Unroll the circuit
#from qiskit.transpiler import PassManager
#from qiskit.transpiler.passes import Unroller
#pass_ = Unroller(['u3', 'cx'])
#pm = PassManager(pass_)
#qc = pm.run(qc) 
#qc.draw('mpl')

#qc.count_ops()



#import qiskit.quantum_info as qi
#unitary = qi.Operator(qc)
#print(unitary)

qc1 = QuantumCircuit(4)
#T1 = np.matmul(unitary_H, unitary)
#V1 = np.matmul(T1, S)

unitary = [[ 1.00000000e+00-6.16343852e-16j,
            2.08166817e-16-7.21644966e-16j,
           -3.53270803e-16-4.93107163e-16j,
           -2.94392336e-17+3.92523115e-17j,
            1.95446967e-32+1.57859836e-16j,
            3.69778549e-32-1.01972331e-32j,
           -1.01972331e-32-1.73472348e-18j,
           -4.10536659e-48-6.84227766e-49j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j],
          [ 3.33066907e-16-6.10622664e-16j,
            7.07106781e-01+7.07106781e-01j,
            5.31833334e-17-1.06089094e-16j,
           -7.13922383e-17+6.01412094e-16j,
           -1.44236986e-18-9.63760726e-19j,
            9.05935082e-17-1.35582766e-16j,
           -1.44236986e-18-9.63760726e-19j,
           -1.74565700e-32+9.30061327e-33j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j],
          [-5.88784672e-16+2.77333912e-32j,
           -3.92523115e-17+3.92523115e-17j,
           -7.07106781e-01-7.07106781e-01j,
            6.66133815e-16-3.88578059e-16j,
            3.56529216e-18-3.79105284e-18j,
           -2.52736856e-18-2.37686144e-18j,
           -1.10524057e-16+1.17522638e-16j,
           -3.69778549e-32-1.84889275e-32j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j],
          [-5.44445200e-17+1.08296884e-17j,
           -4.56771420e-16-4.46829966e-16j,
           -7.77156117e-16+2.22044605e-16j,
            1.38222767e-14+1.00000000e+00j,
            3.90458345e-19+1.69020940e-18j,
           -6.35556882e-33+1.09778007e-32j,
            3.90458345e-19+1.69020940e-18j,
           -1.48738427e-16+3.43603344e-17j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j],
          [-1.66533454e-16+2.03944661e-32j,
           -2.77555756e-17-2.46519033e-32j,
            1.74315280e-32+1.01972331e-32j,
            1.38777878e-17-6.93889390e-18j,
            1.34024056e-15+1.00000000e+00j,
           -7.21644966e-16-2.08166817e-16j,
           -4.98013702e-16+3.53270803e-16j,
            3.92523115e-17-9.81307787e-18j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j],
          [ 2.46519033e-32-3.69778549e-32j,
           -9.60882041e-17-1.30966779e-16j,
           -4.10536659e-48+1.36845553e-48j,
           -1.94060638e-18+1.53940018e-17j,
            6.10622664e-16+3.33066907e-16j,
            7.07106781e-01-7.07106781e-01j,
            2.30963598e-17+3.60861489e-17j,
           -6.34049198e-16-9.31996541e-17j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j],
          [ 8.71576399e-33-5.33686624e-49j,
           -9.50744577e-18+1.01094742e-17j,
            1.21313691e-16+1.14089349e-16j,
            1.48631971e-17+4.45270865e-18j,
           -1.96261557e-17-5.49532361e-16j,
           -1.96261557e-17-7.85046229e-17j,
           -7.07106781e-01+7.07106781e-01j,
            3.33066907e-16+5.55111512e-16j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j],
          [-1.36845553e-48+1.36845553e-48j,
           -2.70433504e-17+6.24733352e-18j,
           -2.46519033e-32+0.00000000e+00j,
            4.64285787e-18-1.64870265e-16j,
           -1.09037079e-17+1.63185521e-17j,
           -3.81481738e-16+4.29623179e-16j,
            2.77555756e-16+8.32667268e-16j,
           -1.00000000e+00+1.48492330e-14j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j],
          [ 0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
           -7.07106781e-01+7.07106781e-01j,
           -7.21644966e-16-8.88178420e-16j,
            3.92523115e-16+1.17756934e-16j,
            3.92523115e-17-9.86076132e-32j,
           -2.49800181e-16+4.16333634e-17j,
            4.62223187e-33-2.77555756e-17j,
           -1.38777878e-17+1.38777878e-17j,
           -1.69953884e-33-1.38777878e-17j],
          [ 0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
           -2.49800181e-16-1.11022302e-15j,
           -1.00000000e+00+1.42941214e-14j,
           -7.39557099e-32+7.39557099e-32j,
           -2.82608171e-16-4.24874511e-16j,
           -2.69591305e-17+1.53678319e-17j,
            1.92281088e-16+1.61859462e-16j,
            3.82887304e-18+1.92490447e-17j,
            1.15389588e-17+7.71008581e-18j],
          [ 0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
           -3.14018492e-16+2.55140025e-16j,
           -7.39557099e-32-4.93038066e-32j,
           -1.00000000e+00+2.35749256e-15j,
           -1.94289029e-16+1.11022302e-15j,
           -2.97263942e-17-8.90541731e-18j,
            1.54074396e-33+1.54074396e-32j,
            1.42134668e-16+1.52721161e-16j,
           -9.50744577e-18+1.01094742e-17j],
          [ 0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
           -1.63185521e-17-1.09037079e-17j,
           -3.91941329e-16+2.13632986e-17j,
            7.21644966e-16-9.43689571e-16j,
            7.07106781e-01+7.07106781e-01j,
            3.01670172e-17+7.27434169e-18j,
            2.70433504e-17-6.24733352e-18j,
            4.36886924e-17+4.15067492e-18j,
            5.93496685e-17+2.56911829e-16j],
          [ 0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            2.22044605e-16-2.52317503e-32j,
            3.08148791e-33+2.08166817e-17j,
            1.74315280e-32-1.38777878e-17j,
           -7.70371978e-33-2.77555756e-17j,
           -7.07106781e-01-7.07106781e-01j,
           -9.43689571e-16+8.32667268e-16j,
            1.17756934e-16-4.31775426e-16j,
            9.86076132e-32+1.96261557e-17j],
          [ 0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
           -1.54201716e-17+2.30779177e-17j,
            1.59552478e-16+1.73372251e-16j,
            1.15389588e-17+7.71008581e-18j,
           -1.54201716e-17+2.30779177e-17j,
            1.33226763e-15-2.49800181e-16j,
            1.49047441e-14+1.00000000e+00j,
           -1.09777275e-17+8.70816242e-17j,
            2.99666919e-16-2.71852502e-16j],
          [ 0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
           -3.92338400e-17+1.20405693e-18j,
           -4.75372289e-18+5.05473712e-18j,
            9.27913535e-17+1.44417772e-16j,
            9.50744577e-18-1.01094742e-17j,
           -3.14018492e-16-2.74766180e-16j,
           -9.86076132e-32+6.16297582e-32j,
            3.06403417e-15+1.00000000e+00j,
           -1.16573418e-15-2.22044605e-16j],
          [ 0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
            0.00000000e+00+0.00000000e+00j,
           -6.24733352e-18-2.70433504e-17j,
            5.25248674e-17-1.92555046e-17j,
           -9.37100028e-18-4.05650256e-17j,
           -1.98116510e-17-2.09072462e-16j,
           -1.08296884e-17-5.44445200e-17j,
           -1.31140573e-16+4.78874914e-16j,
           -9.43689571e-16-6.66133815e-16j,
           -7.07106781e-01+7.07106781e-01j]]

#T1 = np.matmul(unitary_H, unitary)
#V1 = np.matmul(T1, S)

#unitary_d = unitary.diagonal()
#unitary_d = unitary_d.tolist()

qc1.h(0)
qc1.h(1)
qc1.h(2)
qc1.h(3)
qc1.iso(unitary, [0, 1, 2, 3], [])
qc1.h(0)
qc1.h(1)
qc1.h(2)
qc1.h(3)

qc1 = transpile(qc1, basis_gates = ['u3', 'cx'], optimization_level=2)
qc1.count_ops()


#from qiskit import QuantumCircuit, BasicAer, execute
#from qiskit.visualization import plot_state_paulivec
#%matplotlib inline
#backend = BasicAer.get_backend('statevector_simulator')
#job = execute(qc, backend).result()
#plot_state_paulivec(job.get_statevector(qc), color = 'midnightblue')
# obtain gates
#gates=new_circuit.count_ops()


# In[148]:


##### check your quantum circuit by running the next line
check_circuit(qc1)


# You can check whether your circuit is valid before submitting it with `check_circuit(qc)`. Once you have a valid solution, please submit it by running the following cell (delete the `#` before `submit_circuit`). You can re-submit at any time.
# 

# In[140]:


# Send the circuit as the final answer, can re-submit at any time
submit_circuit(qc1) 

